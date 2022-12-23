import numba as nb
import numpy as np

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
    #tr2 = 0
    #tr3 = 0
    for n in nb.prange(n_e):
        l = evals[n]
        tr1 += 1 / l * (evecs[3 * i, n] * evecs[3 * j, n] + evecs[3 * i + 1, n] * evecs[3 * j + 1, n] + evecs[
            3 * i + 2, n] * evecs[3 * j + 2, n])
        # tr2 += 1 / l * (evecs[3 * i, n] * evecs[3 * i, n] + evecs[3 * i + 1, n] * evecs[3 * i + 1, n] + evecs[
        #     3 * i + 2, n] * evecs[3 * i + 2, n])
        # tr3 += 1 / l * (evecs[3 * j, n] * evecs[3 * j, n] + evecs[3 * j + 1, n] * evecs[3 * j + 1, n] + evecs[
        #     3 * j + 2, n] * evecs[3 * j + 2, n])
    cov = tr1#  / np.sqrt(tr2 * tr3)
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
    data = np.exp(-sigma*data ** 2)
    sims = sparse.coo_matrix((data, (d.row, d.col)), shape=d.shape)
    return sims