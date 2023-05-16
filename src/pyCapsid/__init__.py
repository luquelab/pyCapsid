

def read_config(file_path):
    import toml
    params_dict = toml.load(file_path)

    return params_dict



def run_capsid(params_path):

    params_dict = read_config(params_path)


    from .PDB import getCapsid
    pdb = params_dict['PDB']['pdb']
    capsid, calphas, coords, bfactors, chain_starts, title = getCapsid(**params_dict['PDB'])

    from .CG import buildENMPreset

    kirch, hessian = buildENMPreset(coords, **params_dict['CG'])

    from .NMA import modeCalc
    evals, evecs = modeCalc(hessian, **params_dict['NMA'])

    from .NMA import fitCompareBfactors
    evals_scaled, evecs_scaled = fitCompareBfactors(evals, evecs, bfactors, pdb, **params_dict['b_factors'])

    from .NMA import calcDistFlucts
    from .QRC import findQuasiRigidClusters

    dist_flucts = calcDistFlucts(evals_scaled, evecs, coords)


    labels, score, residue_scores = findQuasiRigidClusters(pdb, dist_flucts, **params_dict['QRC'])

    from pyCapsid.VIS import chimeraxViz
    chimeraxViz(labels, pdb, **params_dict['VIS'])