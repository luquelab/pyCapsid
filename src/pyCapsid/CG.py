"""Module with functions for building hessian matrices of different elastic network models."""

import numba as nb
import numpy as np

def buildENMPreset(coords, preset='ANM', **kwargs):
    """Builds a hessian matrix representing an ENM based on one of several presets.

    :param ndarray coords: Cartesian coordinates of alpha carbon atoms(or choice of alternate representation).
    :param str preset: The specific model preset to use. Only accepts the following values:
        - 'ANM': Anisotropic Network Model with a cutoff of 15Å and no distance weighting.
        - 'GNM': Gaussian Network Model with a cutoff of 7.5Å and no distance weighting.
        - 'U-ENM': Unified Elastic Network Model with a cutoff of 7.5Å and and f_anm parameter of 0.1. [ref]
        - 'bbENM': Backbone-enhanced Elastic Network Model with a cutoff of 7.5Å and no distance weighting.
        - 'betaENM': Not yet implemented
    :return: A tuple of sparse matrices. The kirchoff matrix and the hessian matrix
    :rtype: (scipy.sparse.csr_matrix, scipy.sparse.csr_matrix)
    """
    model_presets = ['ANM', 'GNM', 'U-ENM', 'bbENM']
    print('Building hessian for model preset: ', preset)
    if preset == 'ANM':
        cutoff = 15
        return buildENM(coords, cutoff=cutoff, **kwargs)
    elif preset == 'GNM':
        cutoff = 7.5
        gnm=True
        return buildENM(coords, cutoff=cutoff, gnm=gnm, **kwargs)
    elif preset == 'U-ENM':
        cutoff = 7.5
        fanm = 0.1
        return buildENM(coords, cutoff=cutoff, fanm=fanm, **kwargs)
    elif preset == 'bbENM':
        cutoff = 7.5
        l_backbone=1
        k_backbone = 10
        if 'chain_starts' not in kwargs:
            raise ValueError("No chain information provided. Indices of chain starts must be provided as chain_starts")
        chain_starts = kwargs.pop('chain_starts')
        return buildENM(coords, cutoff=cutoff, chain_starts=chain_starts, l_backbone=l_backbone, k_backbone=k_backbone, **kwargs)
    else:
        raise ValueError("Invalid model preset. Expected one of: %s" % model_presets)



def buildENM(coords, cutoff=10, gnm=False, fanm=1, wfunc='power', base_dist=1, d_power=0, backbone=False, k_backbone=1,
             l_backbone=1, chain_starts=None, save_hessian=False, save_kirchhoff=False, save_cg_path='./'):
    """Builds a hessian matrix representing an ENM based on the provided parameters.

        :param coords: Cartesian of alpha carbon (or choice of representation) atoms
        :param float cutoff: Cutoff distance for long range interactions in the ENM
        :param gnm: Whether to use only the kirchhoff matrix (Isotropic GNM), otherwise use full hessian
        :param float fanm: Parameter representing degree of anisotropy for U-ENM
        :param wfunc: Weight function to assign spring constant based on distance. Use 'power' or 'exp'
        :param base_dist: In wfunc, divide distance by the base_distance
        :param d_power: In wfunc, use this power of the distance
        :param backbone: Whether to use stronger interactions for residues connected along the backbone
        :param k_backbone: Relative strength of backbone interaction
        :param l_backbone: How many steps along the backbone to give stronger interactions
        :param chain_starts: Used for defining backbone interactions
        :return: A tuple of sparse matrices. The kirchhoff matrix and the hessian matrix
        :rtype: (scipy.sparse.csr_matrix, scipy.sparse.csr_matrix)
        """
    params = locals()
    params.pop('coords')
    print('Model parameters: ', params)

    import numpy as np
    from scipy import sparse
    from sklearn.neighbors import BallTree, radius_neighbors_graph


    n_atoms = coords.shape[0]
    dof = n_atoms * 3

    print(f'Finding neighbors within {cutoff}Å')
    tree = BallTree(coords)
    distGraph = radius_neighbors_graph(tree, cutoff, mode='distance', n_jobs=-1)
    dists = distGraph.tocoo().copy()
    dists.sum_duplicates()

    print('Building kirchhoff matrix')
    kirch = kirchGamma(dists, gfunc=wfunc, bd=base_dist, d2=d_power).tolil()

    if backbone:
        print('Adding backbone terms')
        try:
            c = chain_starts[0]
        except:
            print('Must provide an array of chain starts')
        buildBackbone(l_backbone, k_backbone, kirch, chain_starts)

    dg = np.array(kirch.sum(axis=0))
    kirch.setdiag(-dg[0])
    kirch = kirch.tocsr()

    if not gnm:
        print('Building hessian matrix')
        kc = kirch.tocoo().copy()
        hData = hessCalc(kc.row, kc.col, kc.data, coords)
        indpt = kirch.indptr
        inds = kirch.indices
        hessian = sparse.bsr_matrix((hData, inds, indpt), shape=(dof, dof)).tocsr()
        hessian = fanm * hessian + (1 - fanm) * sparse.kron(kirch, np.identity(3))
    else:
        hessian = kirch.copy()

    print('Done building model')

    from sklearn.utils.validation import check_symmetric
    check_symmetric(hessian, raise_warning=True, tol=1e-5)
    check_symmetric(kirch, raise_warning=True, tol=1e-5)
    hessian.eliminate_zeros()
    kirch.eliminate_zeros()

    if save_kirchhoff:
        file = save_cg_path + "kirchhoff"
        sparse.save_npz(file, kirch, compressed=True)

    if save_hessian:
        file = save_cg_path + "hessian"
        sparse.save_npz(file, hessian, compressed=True)

    return kirch, hessian


def kirchGamma(dists, **kwargs):
    """ Calculates the kirchhoff matrix (spring constant matrix) from a sparse distance matrix.

    :param dists: Sparse distance matrix
    :param kwargs:
    :return: The sparse kirchhoff matrix.
    """
    if kwargs['gfunc'] == 'power':
        dists.data = -1 / ((dists.data / kwargs['bd']) ** kwargs['d2'])
    elif kwargs['gfunc'] == 'exp':
        dists.data = -np.exp(((dists.data / kwargs['bd']) ** kwargs['d2']))
    else:
        dists.data = -np.ones_like(dists.data)

    return dists


def buildBackbone(bblen, bbgamma, kirch, chain_starts):
    """ Modifies the kirchhoff matrix to include backbone terms

    :param bblen: The number of backbone neighbors in either direction with strengthened spring constants.
    :param bbgamma: The relative strength to use for the backbone interactions.
    :param kirch: The sparse kirchhoff matrix
    :param chain_starts: Indices where protein chains begin/end.
    """
    for i in range(len(chain_starts) - 1):
        start = chain_starts[i]
        stop = chain_starts[i + 1]
        for j in range(stop - start):
            ij = start + j
            for k in range(1, bblen + 1):
                neighbor = ij - k
                if neighbor < start:
                    continue
                kirch[ij, neighbor] = -bbgamma / k
                kirch[neighbor, ij] = -bbgamma / k


@nb.njit()
def hessCalc(row, col, kGamma, coords):
    """
    :param row:
    :param col:
    :param kGamma:
    :param coords:
    :return:
    """
    hData = np.zeros((row.shape[0], 3, 3))
    dinds = np.nonzero(row == col)[0]
    for k, (i, j) in enumerate(zip(row, col)):
        if i == j:
            continue
        dvec = coords[j] - coords[i]
        d2 = np.dot(dvec, dvec)
        g = kGamma[k]
        hblock = np.outer(dvec, dvec) * (g / d2)
        hData[k, :, :] = hblock
        hData[dinds[i]] += -hblock
        # hData[dinds[j]] += -hblock
    return hData


# @nb.njit()
# def cooOperation(row, col, data, func, arg):
#     r = np.copy(data)
#     for n, (i, j, v) in enumerate(zip(row, col, data)):
#         if i == j:
#             continue
#         r[n] = func(i, j, v, arg)
#     return r


# def betaCarbonModel(calphas):
#     """
#
#     :param calphas:
#     :return:
#     """
#     from sklearn.neighbors import BallTree, radius_neighbors_graph
#     from scipy import sparse
#     coords = calphas.getCoords()
#     n_atoms = coords.shape[0]
#     n_asym = int(n_atoms / 60)
#
#     print('# Atoms ', n_atoms)
#     dof = n_atoms * 3
#     betaCoords, bvals, nobetas, aBetas = buildBetas(coords, calphas)
#     print(bvals.shape)
#     print(betaCoords.shape)
#
#     tree = BallTree(coords)
#     kirch = radius_neighbors_graph(tree, cutoff, mode='distance', n_jobs=-1)
#     kc = kirch.tocoo()
#     kc.sum_duplicates()
#     kc.eliminate_zeros()
#     bfactors = calphas.getBetas()
#     kirch = kirchGamma(kc, d2=d2)
#     if backboneConnect:
#         kbb = backboneStrength
#         kirch = backbonePrody(calphas, kirch.tolil(), kbb, s=bblen)
#     print('nonzero values: ', kirch.nnz)
#     dg = np.array(kirch.sum(axis=0))
#     kirch.setdiag(-dg[0])
#     kirch.sum_duplicates()
#     kirch = kirch.tocsr()
#
#     btree = BallTree(betaCoords)
#     betaKirch = radius_neighbors_graph(btree, cutoff, mode='distance', n_jobs=-1).tocoo()
#     betaKirch.eliminate_zeros()
#     betaKirch = kirchGamma(betaKirch.tocoo(), d2=d2)
#     dg = np.array(betaKirch.sum(axis=0))
#     betaKirch.setdiag(-dg[0])
#     betaKirch = betaKirch.tocsr()
#
#     from scipy.spatial import KDTree
#     akdtree = KDTree(coords)
#     bkdtree = KDTree(betaCoords)
#     abKirch = akdtree.sparse_distance_matrix(bkdtree, cutoff).tocoo()
#     abKirch = kirchGamma(abKirch.tocoo(), d2=d2)
#     dg = np.array(abKirch.sum(axis=0))
#     abKirch.setdiag(-dg[0])
#     abKirch.eliminate_zeros()
#     abKirch = abKirch.tocoo()
#
#     print('aa Hess')
#     kc = kirch.tocoo().copy()
#     kc.sum_duplicates()
#     haaData = hessCalc(kc.row, kc.col, kirch.data, coords)
#     indpt = kirch.indptr
#     inds = kirch.indices
#     haa = sparse.bsr_matrix((haaData, inds, indpt), shape=(dof, dof)).tolil()
#     haa = fanm * haa + (1 - fanm) * sparse.kron(kirch, np.identity(3))
#     haa = haa.tocsr()
#     haa.eliminate_zeros()
#     haa.sum_duplicates()
#
#     print('bb Hess')
#     kbc = betaKirch.tocoo().copy()
#     kbc.sum_duplicates()
#     row, col = (kbc.row, kbc.col)
#     hbrow, hbcol, hbdata, kbrow, kbcol, kbdata = bbHess(row, col, betaCoords, nobetas, bvals, aBetas, betaKirch.data,
#                                                         fanm)
#     hbb = sparse.coo_matrix((hbdata, (hbrow, hbcol)), shape=(dof, dof)).tocsr()
#     kbb = sparse.coo_matrix((kbdata, (kbrow, kbcol)), shape=(n_atoms, n_atoms)).tocsr()
#     hbb.eliminate_zeros()
#     hbb.sum_duplicates()
#     kbb.eliminate_zeros()
#     kbb.sum_duplicates()
#
#     print('ab Hess')
#     # abKirch = kirchGamma(abKirch.tocoo(), bfactors, d2=d2).tocoo()
#     abKirch.eliminate_zeros()
#     abKirch.sum_duplicates()
#     row, col, dat = (abKirch.row, abKirch.col, abKirch.data)
#
#     hrow, hcol, habData, krow, kcol, kabData = abHessOnly(row, col, dat, betaCoords, coords, nobetas, bvals, aBetas,
#                                                           fanm)
#     hab = sparse.coo_matrix((habData, (hrow, hcol)), shape=(dof, dof)).tocsr()
#     kab = sparse.coo_matrix((kabData, (krow, kcol)), shape=(n_atoms, n_atoms)).tocsr()
#
#     from sklearn.utils.validation import check_symmetric
#     check_symmetric(haa, raise_exception=True, tol=1e-5)
#     check_symmetric(hbb, raise_exception=True, tol=1e-5)
#     check_symmetric(hab, raise_exception=True, tol=1e-5)
#     print(kirch.shape, kbb.shape, kab.shape)
#     hess = aaGamma * haa + bbGamma * hbb + abGamma * hab
#
#     kirchoff = kirch * aaGamma + kbb * bbGamma + kab * abGamma
#     check_symmetric(hess, raise_warning=True, tol=1e-5)
#     print(hess.data)
#     hess.eliminate_zeros()
#
#     start = 3 * 0
#     stop = 3 * 30
#
#     fig, ax = plt.subplots()
#     mat = ax.matshow(haa[start:stop, start:stop].todense())
#     plt.colorbar(mat, ax=ax)
#     plt.show()
#     checkHessStrength(kirch, haa)
#
#     fig, ax = plt.subplots()
#
#     mat = ax.matshow(hbb[start:stop, start:stop].todense())
#     plt.colorbar(mat, ax=ax)
#     plt.show()
#     checkHessStrength(kbb, hbb)
#
#     fig, ax = plt.subplots()
#     mat = ax.matshow(hab[start:stop, start:stop].todense())
#     plt.colorbar(mat, ax=ax)
#     plt.show()
#     checkHessStrength(kab, hab)
#
#     fig, ax = plt.subplots()
#     mat = ax.matshow(hess[start:stop, start:stop].todense())
#     plt.colorbar(mat, ax=ax)
#     plt.show()
#     checkHessStrength(kirchoff, hess)
#
#     return kirchoff, hess
#
#
# def checkHessStrength(k, h):
#     """
#
#     :param k:
#     :param h:
#     """
#     hd = h.diagonal().reshape((-1, 3))
#     kd = k.diagonal()
#     hd = np.sum(hd, axis=-1)
#     print(np.min(hd))
#     print(hd)
#     print(kd)
#     print('total connections', np.sum(hd), np.sum(kd))
#     print('matrix sum', h.sum())
#
#
# # @nb.njit()
# def buildBetas(coords, calphas):
#     """
#
#     :param coords:
#     :param calphas:
#     :return:
#     """
#     na = calphas.numAtoms()
#     noBetas = []
#     aBetas = []
#     bvals = []
#     bcoords = np.zeros_like(coords)
#     count = 0
#     chid = 0
#     ch0 = calphas[0].getChid()
#     seg0 = calphas[0].getSegname()
#     nch = 0
#     for i, ca in enumerate(calphas.iterAtoms()):
#         if ca.getResname() == 'GLY':
#             # bcoords[i, :] = coords[i, :]
#             noBetas.append(i)
#             bvals.append(0)
#             # print('nobeta', i)
#             count += 1
#         elif count == 0 or i == na - 1:
#             # bcoords[i, :] = coords[i, :]
#             aBetas.append(i)
#             # noBetas.append(i)
#             # print('abeta', i)
#             # b = 1/np.linalg.norm(coords[i, :])
#             bvals.append(1)
#             count += 1
#         elif ca.getChid() == ch0 and ca.getSegname() == seg0 and calphas[i - 1].getChid() == ch0 and calphas[
#             i - 1].getSegname() == seg0 and calphas[i + 1].getChid() == ch0 and calphas[i + 1].getSegname() == seg0:
#             r = 2 * coords[i, :] - coords[i - 1, :] - coords[i + 1, :]
#             b = 1 / np.linalg.norm(r)
#             bvals.append(b)
#             bcoords[i, :] = coords[i, :] + 3.0 * r * b
#             count += 1
#         else:
#             bcoords[i, :] = coords[i, :]
#             aBetas.append(i)
#             bvals.append(1)
#             chid += 1
#             ch0 = ca.getChid()
#             seg0 = ca.getSegname()
#             count = 1
#             nch += 1
#     print(len(aBetas))
#
#     return bcoords, np.array(bvals), np.array(noBetas), np.array(aBetas)
#
#
# @nb.njit()
# def bbHessOnly(row, col, bcoords, nobetas, kGamma):
#     """
#
#     :param row:
#     :param col:
#     :param bcoords:
#     :param nobetas:
#     :param kGamma:
#     :return:
#     """
#     hData = np.zeros((row.shape[0], 3, 3))
#     dinds = np.nonzero(row == col)[0]
#     for k, (i, j) in enumerate(zip(row, col)):
#         if abs(i - j) <= 1 or np.any(nobetas == i) or np.any(nobetas == j):
#             continue
#
#         dvec = bcoords[j] - bcoords[i]
#         d2 = np.dot(dvec, dvec)
#         g = kGamma[k]
#         hblock = np.outer(dvec, dvec) * (g / d2)  # + (1-fanm)*np.identity(3)*g
#         hData[k, :, :] = hblock
#         hData[dinds[i]] += -hblock / 2
#         hData[dinds[j]] += -hblock / 2
#
#     return hData
#
#
# @nb.njit()
# def bbHess(row, col, bcoords, nobetas, bvals, abetas, kGamma, fanm):
#     """
#
#     :param row:
#     :param col:
#     :param bcoords:
#     :param nobetas:
#     :param bvals:
#     :param abetas:
#     :param kGamma:
#     :param fanm:
#     :return:
#     """
#     nrow = []
#     ncol = []
#     hData = []
#     krow = []
#     kcol = []
#     kData = []
#     n = bcoords.shape[0]
#
#     for k, (i, j) in enumerate(zip(row, col)):
#         if np.any(nobetas == j) or np.any(nobetas == i) or abs(i - j) <= 1:
#             continue
#         rvec = bcoords[i] - bcoords[j]
#         d2 = np.dot(rvec, rvec)
#         g = kGamma[k]
#         bi = bvals[i]
#         bj = bvals[j]
#         hblock = -(fanm * np.outer(rvec, rvec) + (1 - fanm) * np.identity(3)) * (
#                 g / d2) / 2  # + (1-fanm)*np.identity(3)*g
#         ijs = (3 * i, 3 * i + 3, 3 * i - 3, 3 * j, 3 * j + 3, 3 * j - 3)
#
#         if np.any(j == abetas):
#             cjs = (-1, 0, 0)
#         else:
#             cjs = (-(1 + 6 * bj), 3 * bj, 3 * bj)
#
#         if np.any(i == abetas):
#             cis = (1, 0, 0)
#         else:
#             cis = ((1 + 6 * bi), -3 * bi, -3 * bi)
#         cs = cis + cjs
#
#         for ci, ii in zip(cs, ijs):
#             if ci == 0:
#                 continue
#             for cj, jj in zip(cs, ijs):
#                 if cj == 0:
#                     continue
#                 krow.append(int(ii / 3))
#                 kcol.append(int(jj / 3))
#                 kData.append(-ci * cj * g)
#                 for l in range(3):
#                     for m in range(3):
#                         nrow.append(ii + l)
#                         ncol.append(jj + m)
#                         hData.append(ci * cj * hblock[l, m])
#
#     return nrow, ncol, hData, krow, kcol, kData
#
#
# @nb.njit()
# def abHessOnly(row, col, kGamma, bcoords, coords, nobetas, bvals, abetas, fanm):
#     """
#
#     :param row:
#     :param col:
#     :param kGamma:
#     :param bcoords:
#     :param coords:
#     :param nobetas:
#     :param bvals:
#     :param abetas:
#     :param fanm:
#     :return:
#     """
#     nrow = []
#     ncol = []
#     hData = []
#
#     krow = []
#     kcol = []
#     kData = []
#     for k, (i, j) in enumerate(zip(row, col)):
#         if np.any(nobetas == j) or abs(i - j) <= 1:
#             continue
#         rvec = coords[i] - bcoords[j]
#         d2 = np.dot(rvec, rvec)
#         g = kGamma[k]
#         b = bvals[j]
#         hblock = -(fanm * np.outer(rvec, rvec) + (1 - fanm) * np.identity(3)) * (g / d2)
#
#         ijs = (3 * i, 3 * j, 3 * j + 3, 3 * j - 3)
#
#         if np.any(j == abetas):
#             cjs = (-1., 0, 0)
#         else:
#             cjs = (-(1 + 6 * b), 3 * b, 3 * b)
#
#         cis = (1.,)
#
#         cs = cis + cjs
#
#         for ci, ii in zip(cs, ijs):
#             if ci == 0:
#                 continue
#             for cj, jj in zip(cs, ijs):
#                 if cj == 0.:
#                     continue
#                 krow.append(int(ii / 3))
#                 kcol.append(int(jj / 3))
#                 kData.append(-ci * cj * g)
#                 for l in range(3):
#                     for m in range(3):
#                         nrow.append(ii + l)
#                         ncol.append(jj + m)
#                         hData.append(ci * cj * hblock[l, m])
#
#     return nrow, ncol, hData, krow, kcol, kData
