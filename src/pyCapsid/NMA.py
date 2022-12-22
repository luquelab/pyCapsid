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
        epredict = cp.random.rand(n_dim, n_modes + 6, dtype = cp.float32)

        lu = spilu(sparse_gpu, fill_factor=50)  # LU decomposition
        M = LinearOperator(mat.shape, lu.solve)
        print('gpu eigen')

        evals, evecs = clobpcg(sparse_gpu, epredict, M=M,  largest=False, tol=0, verbosityLevel=0)
        if model=='anm':
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
        sparse_gpu_shifted = cp.sparse.csr_matrix((mat - sigma*sparse.eye(mat.shape[0])).astype(cp.float32))
        #mat_shifted = sparse.csc_matrix(mat - sigma*sparse.eye(mat.shape[0]))
        #lu = sparse.linalg.splu(mat_shifted)
        A_gpu_LU = cpsp_la.splu(sparse_gpu_shifted)  # LU decomposition
        #A_gpu_LO = cpsp_la.LinearOperator(mat_shifted.shape, lu.solve)  # Linear Operator
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