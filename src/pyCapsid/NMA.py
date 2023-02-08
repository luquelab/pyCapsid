"""Module with functions for calculating the normal modes and frequencies of a given hessian. Eigenvalue/vector functions.
Also contains functions for calculating mechanical properties from NMA results. I.e. squared fluctuations, distance
fluctuations, compressibility, collectivity etc."""
import numba as nb
import numpy as np


def modeCalc(hess, n_modes=200, eigmethod='eigsh', is3d=True):
    """Calculate the 'n_modes' lowest frequency modes of the system by calculating the smallest eigenvalues and eigenvectors
    of the hessian matrix.

    :param hess: Sparse hessian matrix
    :param n_modes: Integer number of low-frequency modes to calculate.
    :param eigmethod: Choice of method for solving the eigenvalue problem.
    :returns:
    """
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

        epredict = np.random.rand(n_dim, n_modes + 6)
        evals, evecs = lobpcg(hess, epredict, largest=False, tol=0, maxiter=n_dim)
        evals = evals[6:]
        evecs = evecs[:, 6:]
    elif eigmethod == 'lobcuda':
        import cupy as cp
        from cupyx.scipy.sparse.linalg import lobpcg as clobpcg
        from cupyx.scipy.sparse.linalg import LinearOperator, spilu

        sparse_gpu = cp.sparse.csr_matrix(hess.astype(cp.float32))
        epredict = cp.random.rand(n_dim, n_modes + 6, dtype=cp.float32)

        lu = spilu(sparse_gpu, fill_factor=50)  # LU decomposition
        M = LinearOperator(hess.shape, lu.solve)
        print('gpu eigen')

        evals, evecs = clobpcg(sparse_gpu, epredict, M=M, largest=False, tol=0, verbosityLevel=0)
        if is3d == 'anm':
            evals = cp.asnumpy(evals[6:])
            evecs = cp.asnumpy(evecs[:, 6:])
        else:
            evals = cp.asnumpy(evals[1:])
            evecs = cp.asnumpy(evecs[:, 1:])

    end = time.time()
    print('NMA time: ', end - start)
    return evals, evecs


def calcCovMat(evals, evecs, n_modes, coords, fluct_cutoff, is3d=True):
    """Calculates a sparse covariance matrix from the low frequency modes.

    :param evals:
    :param evecs:
    :param n_modes:
    :param coords:
    :param fluct_cutoff:
    :param is3d:
    :return:
    """
    from scipy import sparse
    from sklearn.neighbors import BallTree, radius_neighbors_graph

    print('Calculating sparse covariance matrix')

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


def calcDistFlucts(evals, evecs, coords, fluct_cutoff=7.5, is3d=True):
    """Calculates a sparse covariance matrix from the low frequency modes.

    :param evals:
    :param evecs:
    :param n_modes:
    :param coords:
    :param fluct_cutoff:
    :param is3d:
    :return:
    """
    from scipy import sparse

    n_modes = evals.shape[0]
    covariance = calcCovMat(evals, evecs, n_modes, coords, fluct_cutoff, is3d)
    c_diag = covariance.diagonal()
    row, col, c_data = (covariance.row, covariance.col, covariance.data)

    print('Calculating sparse distance fluctuation matrix from covariance matrix')
    fluct_data = distFluctFromCov(c_diag, c_data, row, col)

    dist_flucts = sparse.coo_matrix((fluct_data, (row, col)), shape=covariance.shape)

    return dist_flucts


@nb.njit()
def con_c(evals, evecs, row, col):
    """

    :param evals:
    :param evecs:
    :param row:
    :param col:
    :return:
    """
    data = []
    for k in range(row.shape[0]):
        i, j = (row[k], col[k])
        data.append(cov(evals, evecs, i, j))
    return np.array(data)


@nb.njit()
def gCon_c(evals, evecs, row, col):
    """

    :param evals:
    :param evecs:
    :param row:
    :param col:
    :return:
    """
    data = []
    for k in range(row.shape[0]):
        i, j = (row[k], col[k])
        c = gCov(evals, evecs, i, j)
        data.append(c)
    return np.array(data)


@nb.njit(parallel=True)
def cov(evals, evecs, i, j):
    """

    :param evals:
    :param evecs:
    :param i:
    :param j:
    :return:
    """
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
    """

    :param evals:
    :param evecs:
    :param i:
    :param j:
    :return:
    """
    # Calculates the covariance between two residues in GNM
    n_e = evals.shape[0]
    c = 0
    for n in range(n_e):
        l = evals[n]
        c += 1 / l * (evecs[i, n] * evecs[j, n])
    return c


@nb.njit()
def distFluctFromCov(c_diag, c_data, row, col):
    """

    :param c_diag:
    :param c_data:
    :param row:
    :param col:
    :return:
    """
    d_data = []
    for k in range(row.shape[0]):
        i, j = (row[k], col[k])
        d = np.abs(c_diag[i] + c_diag[j] - 2 * c_data[k])
        d_data.append(d)
    return np.array(d_data)


def fluctPlot(d, title, pdb):
    """

    :param d:
    :param title:
    :param pdb:
    """
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


# This is a stupid way to switch this to be in NMA
def fitCompareBfactors(evals, evecs, bfactors, pdb, is3d=True, fitModes=True, plotModes=False, forceIco=True, icotol=0.002):
    """

    :param evals:
    :param evecs:
    :param bfactors:
    :param pdb:
    :param is3d:
    :param fitModes:
    :param plotModes:
    :param forceIco:
    :param icotol:
    :return:
    """
    from .bfactorfit import fitPlotBfactors
    cc, gamma,n_modes = fitPlotBfactors(evals, evecs, bfactors, pdb, is3d=is3d, fitModes=fitModes, plotModes=plotModes, forceIco=forceIco, icotol=icotol)
    r_evals = evals[:n_modes]*gamma
    r_evecs = evecs[:, :n_modes]
    return r_evals, r_evecs