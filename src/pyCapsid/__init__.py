

def read_config(file_path):
    import toml
    params_dict = toml.load(file_path)

    return params_dict



def run_capsid(params_path):

    params_dict = read_config(params_path)


    from PDB import getCapsid
    pdb = params_dict['PDB']['pdb']
    capsid, calphas, coords, bfactors, chain_starts, title = getCapsid(pdb, **params_dict['PDB'])

    from CG import buildENMPreset
    kirch, hessian = buildENMPreset(coords, preset='U-ENM')

    from pyCapsid.NMA import modeCalc
    evals, evecs = modeCalc(hessian)

    from pyCapsid.NMA import fitCompareBfactors
    evals_scaled, evecs_scaled = fitCompareBfactors(evals, evecs, bfactors, pdb, fitModes=False)

    from pyCapsid.NMA import calcDistFlucts
    from pyCapsid.QRC import findQuasiRigidClusters
    import numpy as np

    dist_flucts = calcDistFlucts(evals_scaled, evecs, coords)

    n_cluster_max = 100
    n_range = np.arange(4, n_cluster_max, 2)
    labels, score, residue_scores = findQuasiRigidClusters(pdb, dist_flucts, n_range)

    from pyCapsid.VIS import chimeraxViz
    chimeraxViz(labels, pdb, chimerax_path='C:\\Program Files\\ChimeraX\\bin')