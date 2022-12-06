import matplotlib.pyplot as plt
import numba as nb
import numpy as np
from settings import *

def buildENM(calphas, coords, bfactors, cbeta=False, gfunc = 'power', backbone=False, k_backbone = 1, l_backbone=1):
    from scipy import sparse
    import numpy as np
    from sklearn.neighbors import BallTree, radius_neighbors_graph, kneighbors_graph
    from settings import cbeta, backboneStrength, backboneConnect, bblen, gfunc, baseDistance

    n_atoms = coords.shape[0]
    dof = n_atoms * 3

    tree = BallTree(coords)
    distGraph = radius_neighbors_graph(tree, cutoff, mode='distance', n_jobs=-1)
    dists = distGraph.tocoo().copy()
    dists.sum_duplicates()

    kirch = kirchGamma(dists, bfactors, d2=d2, flexibilities=flexibilities, gfunc=gfunc, bd= baseDistance)

    if backboneConnect:
        kbb = backboneStrength
        kirch = backbonePrody(calphas, kirch.tolil(), kbb, s=bblen)
    print('nonzero values: ', kirch.nnz)
    dg = np.array(kirch.sum(axis=0))
    kirch.setdiag(-dg[0])
    kirch.sum_duplicates()
    kirch = kirch.tocsr()
    #print(kirch.data)
    #print('kirch: ', kirch.diagonal(k=0))

    if model=='anm':
        kc = kirch.tocoo().copy()
        kc.sum_duplicates()
        hData = hessCalc(kc.row, kc.col, kirch.data, coords)
        indpt = kirch.indptr
        inds = kirch.indices
        hessian = sparse.bsr_matrix((hData, inds, indpt), shape=(dof,dof)).tocsr()
        hessian = fanm*hessian + (1-fanm)*sparse.kron(kirch, np.identity(3))
    else:
        hessian = kirch.copy()
    print('done constructing matrix')
    # fig, ax = plt.subplots()
    # mat = ax.matshow(hessian[:int(3 * n_asym / 10), :int(3 * n_asym / 10)].todense())
    # plt.colorbar(mat, ax=ax)
    # plt.show()
    from sklearn.utils.validation import check_symmetric
    check_symmetric(hessian, raise_warning=True, tol=1e-5)
    check_symmetric(kirch, raise_warning=True, tol=1e-5)
    return kirch, hessian

def kirchGamma(dists, **kwargs):

    if kwargs['gfunc'] == 'power':
        dists.data = -1/((dists.data/kwargs['bd'])**kwargs['d2'])
        print(dists.data)
    elif kwargs['gfunc'] == 'exp':
        dists.data = -np.exp(((dists.data/kwargs['bd']) ** kwargs['d2']))
    else:
        dists.data = -np.ones_like(dists.data)



    return dists


# need to test this version
def buildBackbone(bblen, bbgamma, kirch, chain_starts):
    for i in range(len(chain_starts)):
        start = chain_starts[i]
        stop = chain_starts[i+1] - 1
        for j in range(start, stop-bblen):
            ij  = start + j
            for k in range(1,bblen):
                kirch[ij, ij + k] = bbgamma/k
                kirch[ij + k, ij] = bbgamma/k

    return kirch






def buildMassesCoords(atoms):
    print(atoms[0])
    bfs = []
    masses = []

    #print('segments:', atoms.numSegments())

    for chain in atoms.iterChains():
        for res in chain.iterResidues():
            mass = np.sum(res.getMasses())
            masses.append(mass)

            bfactor = np.mean(res.getBetas())
            bfs.append(bfactor)

            #coord = res['CA'].getCoords()
            #coords.append(coord)

    return np.asarray(bfs), np.asarray(masses)



@nb.njit()
def cooOperation(row, col, data, func, arg):
    r = np.copy(data)
    for n, (i, j, v) in enumerate(zip(row, col, data)):
        if i==j:
            continue
        r[n] = func(i,j,v, arg)
    return r



@nb.njit()
def hessCalc(row, col, kGamma, coords):
    hData = np.zeros((row.shape[0], 3, 3))
    dinds = np.nonzero(row==col)[0]
    for k, (i,j) in enumerate(zip(row,col)):
        if i==j:
            continue
        dvec = coords[j] - coords[i]
        d2 = np.dot(dvec, dvec)
        g = kGamma[k]
        hblock = np.outer(dvec, dvec) * (g / d2) # + (1-fanm)*np.identity(3)*g
        hData[k,:,:] = hblock
        hData[dinds[i]] += -hblock
        #hData[dinds[j]] += -hblock
    return hData

def backbonePrody(calphas, kirch, k, s):
    k = -k
    nr = calphas.numAtoms()
    count = 0
    chid = 0
    ch0 = calphas[0].getChid()
    seg0 = calphas[0].getSegname()
    # kirch[0,1:3] = -k
    # kirch[1:3, 0] = -k
    nch = 0
    #kbbs = np.array([k,k/2,k/4])[:s]
    for i, ca in enumerate(calphas.iterAtoms()):
        if count == 0:
            for j in range(2*s):
                kirch[i,i+j+1] = k/(j+1)
                kirch[i+j+1, i] = k / (j+1)
            count += 1
            continue
        elif count<s:
            for j in range(count):
                kirch[i, i-j-1] = k/(j+1)
                kirch[i-j-1, i] = k / (j+1)
            for j in range(2*s - count):
                kirch[i,i+j+1] = k/(j+1)
                kirch[i+j+1, i] = k / (j+1)
            count += 1
            continue
        if ca.getChid() == ch0 and ca.getSegname() == seg0:
            for j in range(s):
                kirch[i,i-j-1] = k/(j+1)
                kirch[i - j - 1,i] = k / (j + 1)
            count += 1
        else:
            chid += 1
            #print(ch0, seg0, 'done')
            ch0 = ca.getChid()
            seg0 = ca.getSegname()
            count = 0
            nch += 1
    print(nch)
    return kirch.tocoo()




def betaCarbonModel(calphas):
    from settings import cutoff, d2, fanm, backboneStrength, backboneConnect, bblen, model
    from sklearn.neighbors import BallTree, radius_neighbors_graph
    from scipy import sparse
    coords = calphas.getCoords()
    n_atoms = coords.shape[0]
    n_asym = int(n_atoms / 60)

    print('# Atoms ', n_atoms)
    dof = n_atoms * 3
    betaCoords, bvals, nobetas, aBetas = buildBetas(coords, calphas)
    print(bvals.shape)
    print(betaCoords.shape)

    tree = BallTree(coords)
    kirch = radius_neighbors_graph(tree, cutoff, mode='distance', n_jobs=-1)
    kc = kirch.tocoo()
    kc.sum_duplicates()
    kc.eliminate_zeros()
    bfactors = calphas.getBetas()
    kirch = kirchGamma(kc, bfactors, d2=d2)
    if backboneConnect:
        kbb = backboneStrength
        kirch = backbonePrody(calphas, kirch.tolil(), kbb, s=bblen)
    print('nonzero values: ', kirch.nnz)
    dg = np.array(kirch.sum(axis=0))
    kirch.setdiag(-dg[0])
    kirch.sum_duplicates()
    kirch = kirch.tocsr()


    btree = BallTree(betaCoords)
    betaKirch = radius_neighbors_graph(btree, cutoff, mode='distance',  n_jobs=-1).tocoo()
    betaKirch.eliminate_zeros()
    betaKirch = kirchGamma(betaKirch.tocoo(), bfactors, d2=d2)
    dg = np.array(betaKirch.sum(axis=0))
    betaKirch.setdiag(-dg[0])
    betaKirch = betaKirch.tocsr()

    from scipy.spatial import KDTree
    akdtree = KDTree(coords)
    bkdtree = KDTree(betaCoords)
    abKirch = akdtree.sparse_distance_matrix(bkdtree, cutoff).tocoo()
    abKirch = kirchGamma(abKirch.tocoo(), bfactors, d2=d2)
    dg = np.array(abKirch.sum(axis=0))
    abKirch.setdiag(-dg[0])
    abKirch.eliminate_zeros()
    abKirch = abKirch.tocoo()

    print('aa Hess')
    kc = kirch.tocoo().copy()
    kc.sum_duplicates()
    haaData = hessCalc(kc.row, kc.col, kirch.data, coords)
    indpt = kirch.indptr
    inds = kirch.indices
    haa = sparse.bsr_matrix((haaData, inds, indpt), shape=(dof, dof)).tolil()
    haa = fanm * haa + (1 - fanm) * sparse.kron(kirch, np.identity(3))
    haa = haa.tocsr()
    haa.eliminate_zeros()
    haa.sum_duplicates()

    print('bb Hess')
    kbc = betaKirch.tocoo().copy()
    kbc.sum_duplicates()
    row, col = (kbc.row, kbc.col)
    hbrow, hbcol, hbdata, kbrow, kbcol, kbdata = bbHess(row, col, betaCoords, nobetas, bvals, aBetas, betaKirch.data, fanm)
    hbb = sparse.coo_matrix((hbdata, (hbrow, hbcol)), shape=(dof, dof)).tocsr()
    kbb = sparse.coo_matrix((kbdata, (kbrow, kbcol)), shape=(n_atoms, n_atoms)).tocsr()
    hbb.eliminate_zeros()
    hbb.sum_duplicates()
    kbb.eliminate_zeros()
    kbb.sum_duplicates()

    print('ab Hess')
    #abKirch = kirchGamma(abKirch.tocoo(), bfactors, d2=d2).tocoo()
    abKirch.eliminate_zeros()
    abKirch.sum_duplicates()
    row, col, dat = (abKirch.row, abKirch.col, abKirch.data)

    hrow, hcol, habData, krow, kcol, kabData = abHessOnly(row, col, dat, betaCoords, coords, nobetas, bvals, aBetas, fanm)
    hab = sparse.coo_matrix((habData, (hrow, hcol)), shape=(dof,dof)).tocsr()
    kab = sparse.coo_matrix((kabData, (krow, kcol)), shape=(n_atoms, n_atoms)).tocsr()


    from sklearn.utils.validation import check_symmetric
    check_symmetric(haa, raise_exception=True, tol=1e-5)
    check_symmetric(hbb, raise_exception=True, tol=1e-5)
    check_symmetric(hab, raise_exception=True, tol=1e-5)
    from settings import aaGamma, bbGamma, abGamma
    print(kirch.shape, kbb.shape, kab.shape)
    hess = aaGamma*haa + bbGamma*hbb + abGamma*hab

    kirchoff = kirch * aaGamma + kbb * bbGamma + kab * abGamma
    check_symmetric(hess, raise_warning=True, tol=1e-5)
    print(hess.data)
    hess.eliminate_zeros()

    start = 3*0
    stop = 3*30

    fig, ax = plt.subplots()
    mat = ax.matshow(haa[start:stop, start:stop].todense())
    plt.colorbar(mat, ax=ax)
    plt.show()
    checkHessStrength(kirch, haa)

    fig, ax = plt.subplots()

    mat = ax.matshow(hbb[start:stop, start:stop].todense())
    plt.colorbar(mat, ax=ax)
    plt.show()
    checkHessStrength(kbb, hbb)

    fig, ax = plt.subplots()
    mat = ax.matshow(hab[start:stop, start:stop].todense())
    plt.colorbar(mat, ax=ax)
    plt.show()
    checkHessStrength(kab, hab)

    fig, ax = plt.subplots()
    mat = ax.matshow(hess[start:stop, start:stop].todense())
    plt.colorbar(mat, ax=ax)
    plt.show()
    checkHessStrength(kirchoff, hess)


    return kirchoff, hess

def checkHessStrength(k, h):
    hd = h.diagonal().reshape((-1,3))
    kd = k.diagonal()
    hd = np.sum(hd,axis=-1)
    print(np.min(hd))
    print(hd)
    print(kd)
    print('total connections', np.sum(hd),np.sum(kd))
    print('matrix sum', h.sum())



# @nb.njit()
def buildBetas(coords, calphas):
    na = calphas.numAtoms()
    noBetas = []
    aBetas = []
    bvals = []
    bcoords = np.zeros_like(coords)
    count = 0
    chid = 0
    ch0 = calphas[0].getChid()
    seg0 = calphas[0].getSegname()
    nch = 0
    for i, ca in enumerate(calphas.iterAtoms()):
        if ca.getResname()=='GLY':
            #bcoords[i, :] = coords[i, :]
            noBetas.append(i)
            bvals.append(0)
            #print('nobeta', i)
            count += 1
        elif count==0 or i==na-1:
            #bcoords[i, :] = coords[i, :]
            aBetas.append(i)
            #noBetas.append(i)
           # print('abeta', i)
            # b = 1/np.linalg.norm(coords[i, :])
            bvals.append(1)
            count += 1
        elif ca.getChid() == ch0 and ca.getSegname() == seg0 and calphas[i-1].getChid() == ch0 and calphas[i-1].getSegname() == seg0 and calphas[i+1].getChid() == ch0 and calphas[i+1].getSegname() == seg0:
            r = 2 * coords[i, :] - coords[i - 1, :] - coords[i + 1, :]
            b = 1/np.linalg.norm(r)
            bvals.append(b)
            bcoords[i, :] = coords[i, :] + 3.0 * r * b
            count += 1
        else:
            bcoords[i, :] = coords[i, :]
            aBetas.append(i)
            bvals.append(1)
            chid += 1
            ch0 = ca.getChid()
            seg0 = ca.getSegname()
            count = 1
            nch += 1
    print(len(aBetas))

    return bcoords, np.array(bvals), np.array(noBetas), np.array(aBetas)




@nb.njit()
def bbHessOnly(row, col, bcoords, nobetas, kGamma):

    hData = np.zeros((row.shape[0], 3, 3))
    dinds = np.nonzero(row == col)[0]
    for k, (i, j) in enumerate(zip(row, col)):
        if abs(i - j) <= 1 or np.any(nobetas==i) or np.any(nobetas==j):
            continue

        dvec = bcoords[j] - bcoords[i]
        d2 = np.dot(dvec, dvec)
        g = kGamma[k]
        hblock = np.outer(dvec, dvec) * (g / d2)  # + (1-fanm)*np.identity(3)*g
        hData[k, :, :] = hblock
        hData[dinds[i]] += -hblock / 2
        hData[dinds[j]] += -hblock / 2

    return hData

@nb.njit()
def bbHess(row, col, bcoords, nobetas, bvals, abetas, kGamma, fanm):

    nrow = []
    ncol = []
    hData = []
    krow = []
    kcol = []
    kData = []
    n = bcoords.shape[0]

    for k, (i,j) in enumerate(zip(row, col)):
        if np.any(nobetas == j) or np.any(nobetas==i) or abs(i-j) <= 1:
            continue
        rvec = bcoords[i] - bcoords[j]
        d2 = np.dot(rvec, rvec)
        g = kGamma[k]
        bi = bvals[i]
        bj = bvals[j]
        hblock = -(fanm*np.outer(rvec, rvec)+ (1-fanm)*np.identity(3)) * (g / d2)/2  # + (1-fanm)*np.identity(3)*g
        ijs = (3*i, 3*i+3, 3*i-3, 3*j, 3*j+3, 3*j-3)



        if np.any(j == abetas):
            cjs = (-1,0,0)
        else:
            cjs = (-(1 + 6 * bj),3 * bj, 3 * bj)

        if np.any(i == abetas):
            cis = (1,0,0)
        else:
            cis = ((1 + 6 * bi),-3 * bi, -3 * bi)
        cs = cis + cjs


        for ci, ii in zip(cs,ijs):
            if ci==0:
                continue
            for cj, jj in zip(cs,ijs):
                if cj==0:
                    continue
                krow.append(int(ii/3))
                kcol.append(int(jj / 3))
                kData.append(-ci*cj*g)
                for l in range(3):
                    for m in range(3):
                        nrow.append(ii+l)
                        ncol.append(jj+m)
                        hData.append(ci*cj*hblock[l,m])

    return nrow, ncol, hData, krow, kcol, kData

@nb.njit()
def abHessOnly(row, col, kGamma, bcoords, coords, nobetas, bvals, abetas, fanm):

    nrow = []
    ncol = []
    hData = []

    krow = []
    kcol = []
    kData = []
    for k, (i,j) in enumerate(zip(row, col)):
        if np.any(nobetas == j) or abs(i-j) <= 1:
            continue
        rvec = coords[i] - bcoords[j]
        d2 = np.dot(rvec, rvec)
        g = kGamma[k]
        b = bvals[j]
        hblock = -(fanm*np.outer(rvec, rvec) + (1-fanm)*np.identity(3)) * (g / d2)

        ijs = (3*i, 3*j, 3*j+3, 3*j-3)


        if np.any(j == abetas):
            cjs = (-1.,0,0)
        else:
            cjs = (-(1 + 6 * b),3 * b, 3 * b)

        cis = (1.,)

        cs = cis + cjs

        for ci, ii in zip(cs,ijs):
            if ci==0:
                continue
            for cj, jj in zip(cs,ijs):
                if cj==0.:
                    continue
                krow.append(int(ii / 3))
                kcol.append(int(jj / 3))
                kData.append(-ci * cj * g)
                for l in range(3):
                    for m in range(3):
                        nrow.append(ii+l)
                        ncol.append(jj+m)
                        hData.append(ci*cj*hblock[l,m])

    return nrow, ncol, hData, krow, kcol, kData

def buildBetaTransform(natoms, bvals, nobetas, aBetas):
    from scipy import sparse
    dof = 3 * natoms
    bT = sparse.lil_matrix((dof, dof))
    for i in range(natoms):
        if i not in nobetas:
            if i in aBetas:
                for j in range(3):
                    ind = 3 * i
                    bT[ind + j, ind + j] = 1
            else:
                b = bvals[i]
                for j in range(3):
                    ind = 3*i
                    bT[ind+j,ind+j] = 1 + 6*b
                    bT[ind + j, ind + j-3] = -3*b
                    bT[ind + j, ind + j + 3] = -3*b
    return bT.tocsr()


def addSulfideBonds(sulfs, kirch):
    sCoords = sulfs.getCoords()
    from sklearn.neighbors import BallTree, radius_neighbors_graph
    tree = BallTree(sCoords)
    adjacency = radius_neighbors_graph(tree, 3.0, n_jobs=-1)
    print(adjacency.nnz)
    sNodes = sulfs.getData('nodeid')
    kirch = kirch.tocsr()
    for i, n in enumerate(sNodes):
        neighbors = np.nonzero(adjacency[i,:]==1)
        print(neighbors)
        kirch[i, neighbors] = kirch[i, neighbors]*100

    return kirch.tocoo()


def addIonBonds(atoms, kirch, dists):
    anion = atoms.select('resname ASP GLU')
    cation = atoms.select('resname ASP GLU')

    for i, at in enumerate(atoms.iterAtoms()):
        if at.getData('resname') == 'ASP or GLU':
            neighbors = dists[i, :] <= 4
            print(np.count_nonzero(neighbors))
            kirch[i, neighbors] = kirch[i, neighbors] * 100

    return kirch

def buildChemENM(capsid):
    from scipy import sparse
    import numpy as np
    from sklearn.neighbors import BallTree, radius_neighbors_graph, kneighbors_graph
    import biotite.structure as struc
    calphas = capsid[capsid.atom_name == 'CA']
    print("Number of protein residues:", struc.get_residue_count(capsid))
    coords = calphas.coord #struc.apply_residue_wise(capsid, capsid.coord, np.average, axis=0)
    n_atoms = coords.shape[0]

    from bondStrengths import detect_disulfide_bonds
    sulfs = detect_disulfide_bonds(capsid)
    print(sulfs.shape[0])

    from bondStrengths import detect_salt_bridges
    salts = detect_salt_bridges(capsid)
    print(salts.shape)

    from bondStrengths import hbondSites
    hbonds = hbondSites(capsid)
    print(hbonds.shape[0])

    print('# Atoms ',n_atoms)
    dof = n_atoms * 3

    tree = BallTree(coords)
    kirch = radius_neighbors_graph(tree, 7.5, mode='distance', n_jobs=-1)
    kc = kirch.tocoo().copy()
    kc.sum_duplicates()
    kirch = kirchChem(kc.tocsr(), hbonds, sulfs, salts)
    # kirch = betaTestKirch(coords)
    print('nonzero values: ', kirch.nnz)
    dg = np.array(kirch.sum(axis=0))
    kirch.setdiag(-dg[0])
    kirch.sum_duplicates()
    print(kirch.data)

    if model=='anm':
        kc = kirch.tocoo().copy()
        kc.sum_duplicates()
        hData = hessCalc(kc.row, kc.col, kirch.data, coords)
        indpt = kirch.indptr
        inds = kirch.indices
        hessian = sparse.bsr_matrix((hData, inds, indpt), shape=(dof,dof)).tocsr()
        hessian = fanm*hessian + (1-fanm)*sparse.kron(kirch, np.identity(3))
    else:
        hessian = kirch.copy()
    print('done constructing matrix')
    return kirch, hessian

def kirchChem(kirch, hbonds, sulfbonds, salts):
    kg = kirch.copy()
    kg.data = -1/kg.data**2

    for ij in hbonds:
        i, j = (ij[0], ij[1])
        kg[i,j] = -10

    for ij in salts:
        i, j = (ij[0], ij[1])
        kg[i,j] = -10

    for ij in sulfbonds:
        i, j = (ij[0], ij[1])
        kg[i,j] = -100

    kg.data = np.where(kg.data<=-1/(3**2), -100, kg.data)
    return kg

def bondAngles(ivec, jvec, kvec):
    rij = jvec - ivec
    rjk = kvec - jvec
    dij = np.linalg.norm(rij)
    djk = np.linalg.norm(rij)
    G = np.dot(rij, rjk) / (dij * djk)
    rblock = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            rblock[i,j]
