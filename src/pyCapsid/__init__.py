

def read_config(file_path):
    import toml
    params_dict = toml.load(file_path)

    return params_dict

def create_directories(params):
    from pathlib import Path
    # pdb = params['PDB']['pdb']
    # results_dir = f'./{pdb}/'
    # Path(results_dir).mkdir(parents=True, exist_ok=True)

    save_suffix = '_path'
    save_params = {}
    for k, v in params.items():
        s = {key: val for key, val in v.items() if key.endswith(save_suffix)}
        save_params.update(s)

    save_params.pop('chimerax_path')
    print(save_params)
    for key, val in save_params.items():

        Path(val).mkdir(parents=True, exist_ok=True)



def run_capsid(params_path):
    from timeit import default_timer as timer
    pycap_start = timer()

    params_dict = read_config(params_path)
    create_directories(params_dict)

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