import numba as nb
import numpy as np

def modeCalc(hess, n_modes, eigmethod, model):
    import time
    print('Calculating Normal Modes')
    start = time.time()

    n_dim = hess.shape[0]

    # if useMass:
    #     from scipy.sparse import diags
    #     if model == 'anm':
    #         masses = np.repeat(masses, 3)
    #
    #     mass = np.tile(masses, 60)
    #     print('mass variance: ', np.std(mass))
    #     print('mass mean: ', np.mean(mass))
    #     print('mass min: ', np.min(mass))
    #     print('mass max: ', np.max(mass))
    #     Mass = diags(mass)
    # else:

    from scipy.sparse import identity
    Mass = identity(n_dim)

    if eigmethod == 'eigsh':
        from scipy.sparse.linalg import eigsh
        evals, evecs = eigsh(hess, k=n_modes, M=Mass, sigma=1e-10, which='LA')
    elif eigmethod == 'lobpcg':
        from scipy.sparse.linalg import lobpcg
        from pyamg import smoothed_aggregation_solver
        diag_shift = 1e-5 * sparse.eye(n_dim)
        mat += diag_shift
        ml = smoothed_aggregation_solver(mat)
        mat -= diag_shift
        M = ml.aspreconditioner()
        epredict = np.random.rand(n_dim, n_modes + 6)
        evals, evecs = lobpcg(mat, epredict, M=M, largest=False, tol=0, maxiter=n_dim)
        evals = evals[6:]
        evecs = evecs[:, 6:]
    elif eigmethod == 'lobcuda':
        import cupy as cp
        from cupyx.scipy.sparse.linalg import lobpcg as clobpcg
        from cupyx.scipy.sparse.linalg import LinearOperator, spilu

        sparse_gpu = cp.sparse.csr_matrix(mat.astype(cp.float32))
        epredict = cp.random.rand(n_dim, n_modes + 6, dtype=cp.float32)

        lu = spilu(sparse_gpu, fill_factor=50)  # LU decomposition
        M = LinearOperator(mat.shape, lu.solve)
        print('gpu eigen')

        evals, evecs = clobpcg(sparse_gpu, epredict, M=M, largest=False, tol=0, verbosityLevel=0)
        if model == 'anm':
            evals = cp.asnumpy(evals[6:])
            evecs = cp.asnumpy(evecs[:, 6:])
        else:
            evals = cp.asnumpy(evals[1:])
            evecs = cp.asnumpy(evecs[:, 1:])
    elif eigmethod == 'eigshcuda':
        import cupy as cp
        import cupyx.scipy.sparse as cpsp
        import cupyx.scipy.sparse.linalg as cpsp_la
        sigma = 1e-10
        sparse_gpu_shifted = cp.sparse.csr_matrix((mat - sigma * sparse.eye(mat.shape[0])).astype(cp.float32))
        # mat_shifted = sparse.csc_matrix(mat - sigma*sparse.eye(mat.shape[0]))
        # lu = sparse.linalg.splu(mat_shifted)
        A_gpu_LU = cpsp_la.splu(sparse_gpu_shifted)  # LU decomposition
        # A_gpu_LO = cpsp_la.LinearOperator(mat_shifted.shape, lu.solve)  # Linear Operator
        A_gpu_LO = cpsp_la.LinearOperator(sparse_gpu_shifted.shape, A_gpu_LU.solve)

        eigenvalues_gpu, eigenstates_gpu = cpsp_la.eigsh(A_gpu_LO, k=n_modes, which='LA', tol=0)

        eigenvalues_gpu = eigenvalues_gpu.get()
        eigenstates_gpu = eigenstates_gpu.get()
        eigenvalues_gpu = (1 + eigenvalues_gpu * sigma) / eigenvalues_gpu
        idx = np.argsort(eigenvalues_gpu)
        eigenvalues_gpu = eigenvalues_gpu[idx]
        evals = cp.asnumpy(eigenvalues_gpu)
        evecs = cp.asnumpy(eigenstates_gpu)

    end = time.time()
    print('NMA time: ', end - start)
    return evals, evecs




def calcCovMat(evals, evecs, n_modes, coords, fluct_cutoff, is3d=True):
    from scipy import sparse
    from sklearn.neighbors import BallTree, radius_neighbors_graph

    n_modes = int(n_modes)
    print(n_modes)
    print('Direct Calculation Method')

    tree = BallTree(coords)
    adj = radius_neighbors_graph(tree, fluct_cutoff, mode='connectivity', n_jobs=-1).tocoo()
    adj.setdiag(1)
    row, col = adj.row, adj.col

    if is3d:
        data = con_c(evals, evecs, row, col)
        covariance = sparse.coo_matrix((data, (row, col)), shape=adj.shape)
    else:
        data = gCon_c(evals, evecs, row, col)
        covariance = sparse.coo_matrix((data, (row, col)), shape=adj.shape)

    return covariance


def calcDistFlucts(evals, evecs, n_modes, coords, fluct_cutoff, is3d=True):
    from scipy import sparse

    covariance = calcCovMat(evals, evecs, n_modes, coords, fluct_cutoff, is3d)
    c_diag = covariance.diagonal()
    row, col, c_data = (covariance.row, covariance.col, covariance.data)
    print('diag: ', c_diag)
    print('data: ', c_data)

    fluct_data = distFluctFromCov(c_diag, c_data, row, col)
    print(fluct_data)

    dist_flucts = sparse.coo_matrix((fluct_data, (row, col)), shape=covariance.shape)

    return dist_flucts


@nb.njit()
def con_c(evals, evecs, row, col):
    data = []
    for k in range(row.shape[0]):
        i, j = (row[k], col[k])
        data.append(cov(evals, evecs, i, j))
    return np.array(data)


@nb.njit()
def gCon_c(evals, evecs, row, col):
    data = []
    for k in range(row.shape[0]):
        i, j = (row[k], col[k])
        c = gCov(evals, evecs, i, j)
        data.append(c)
    return np.array(data)


@nb.njit(parallel=True)
def cov(evals, evecs, i, j):
    # Calculates the covariance between two residues in ANM. Takes the trace of block so no anisotropy info.
    # Commented section would compute normalized covariances
    n_e = evals.shape[0]
    # n_d = evecs.shape[1]
    tr1 = 0
    # tr2 = 0
    # tr3 = 0
    for n in nb.prange(n_e):
        l = evals[n]
        tr1 += 1 / l * (evecs[3 * i, n] * evecs[3 * j, n] + evecs[3 * i + 1, n] * evecs[3 * j + 1, n] + evecs[
            3 * i + 2, n] * evecs[3 * j + 2, n])
        # tr2 += 1 / l * (evecs[3 * i, n] * evecs[3 * i, n] + evecs[3 * i + 1, n] * evecs[3 * i + 1, n] + evecs[
        #     3 * i + 2, n] * evecs[3 * i + 2, n])
        # tr3 += 1 / l * (evecs[3 * j, n] * evecs[3 * j, n] + evecs[3 * j + 1, n] * evecs[3 * j + 1, n] + evecs[
        #     3 * j + 2, n] * evecs[3 * j + 2, n])
    cov = tr1  # / np.sqrt(tr2 * tr3)
    return cov


@nb.njit(parallel=False)
def gCov(evals, evecs, i, j):
    # Calculates the covariance between two residues in GNM
    n_e = evals.shape[0]
    c = 0
    for n in range(n_e):
        l = evals[n]
        c += 1 / l * (evecs[i, n] * evecs[j, n])
    return c


@nb.njit()
def distFluctFromCov(c_diag, c_data, row, col):
    d_data = []
    for k in range(row.shape[0]):
        i, j = (row[k], col[k])
        d = np.abs(c_diag[i] + c_diag[j] - 2 * c_data[k])
        d_data.append(d)
    return np.array(d_data)


def fluctPlot(d, title, pdb):
    import matplotlib.pyplot as plt
    print('Plotting Fluctuation Histogram')
    import matplotlib
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 11}
    matplotlib.rc('font', **font)
    print('Plotting')
    ax.set_ylabel('Density', fontsize=12)
    ax.set_xlabel('$A^{2}$', fontsize=12)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax.legend()
    fig.suptitle(
        'Histogram Of Pairwise Fluctuations: ' + title + ' (' + pdb + ')', fontsize=12)

    ax.hist(d.data, bins='fd', histtype='stepfilled', density=True)
    fig.tight_layout()
    plt.show()


def fluctToSims(d):
    from scipy import sparse
    d_bar = np.mean(np.sqrt(d.data))
    print('RMS distance fluctuations: ', d_bar)
    sigma = 1 / (2 * d_bar ** 2)
    data = d.data
    data = np.exp(-sigma * data ** 2)
    sims = sparse.coo_matrix((data, (d.row, d.col)), shape=d.shape)
    return sims