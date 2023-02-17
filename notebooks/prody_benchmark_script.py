#!/usr/dataB/luquelab/members/ctbrown/miniconda3/envs/pycapsid/bin/python3
#PBS -l nodes=1:ppn=24
#PBS -l walltime=96:00:00

if __name__ == '__main__':
    from timeit import default_timer as timer
    import numpy as np

    from pyCapsid.PDB import getCapsid
    from pyCapsid.ENM import buildENM
    from pyCapsid.NMA import modeCalc
    from prody import ANM, parsePDB

    pdb_files =  ['7kq5', '2e0z']
    cutoff = 10
    n_modes = 200

    for pdb in pdb_files:


        pycap_start = timer()
        capsid, calphas, coords, bfactors, chain_starts, title = getCapsid(pdb, save=False)
        kirch, hessian = buildENM(coords, cutoff=cutoff)
        evals, evecs = modeCalc(hessian, n_modes=n_modes)
        pycap_time = timer() - pycap_start
        print(f'pyCapsid time for {pdb}: {pycap_time}')

        prody_start = timer()
        capsid = parsePDB(pdb, biomol=True)
        calphas = capsid.select('protein and name CA')
        anm = ANM('T. maritima ANM analysis')
        anm.buildHessian(calphas, cutoff=cutoff, sparse=True)
        anm.calcModes(n_modes=n_modes)
        evals_prody = anm.getEigvals()
        evecs_prody = anm.getEigvecs()
        prody_time = timer() - prody_start
        print(f'ProDy time for {pdb}: {prody_time}')
