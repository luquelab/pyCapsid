

def read_config(file_path):
    import toml
    params_dict = toml.load(file_path)

    return params_dict

def create_directories(params_dict):
    from pathlib import Path


    save_suffix = '_path'
    params = params_dict.copy()

    top_path = params['PDB']['save_all_path']

    # This is remarkably hacky but it does its job of creating all the folders for saving results

    for k, v in params.items():
        for key, val in v.items():
            if key.endswith(save_suffix) and not key == 'save_all_path' and not key == 'chimerax_path':
                params_dict[k][key] = top_path + val
                Path(params_dict[k][key]).mkdir(parents=True, exist_ok=True)






def run_capsid(params_path):
    from timeit import default_timer as timer
    pycap_start = timer()

    params_dict = read_config(params_path)
    create_directories(params_dict)
    params_dict['PDB'].pop('save_all_path')

    if params_dict['plotting']['suppress_plots']:
        import matplotlib
        matplotlib.use('Agg')


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

    pycap_time = timer() - pycap_start
    print(f'pyCapsid total execution time for {pdb}: {pycap_time}')

    from pyCapsid.VIS import chimeraxViz
    chimeraxViz(labels, pdb, **params_dict['VIS'])