

def read_config(file_path):
    import toml
    params_dict = toml.load(file_path)

    return params_dict

def create_directories(params_dict):
    from pathlib import Path


    save_suffix = '_path'
    params = params_dict.copy()

    top_path = params['PDB']['save_all_path']
    Path(top_path).mkdir(parents=True, exist_ok=True)

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

    if 'suppress_plots' in params_dict['plotting'].keys():
        if params_dict['plotting']['suppress_plots']:
            import matplotlib
            matplotlib.use('Agg')


    from .PDB import getCapsid
    pdb = params_dict['PDB']['pdb']
    capsid, calphas, asymmetric_unit, coords, bfactors, chain_starts, title = getCapsid(**params_dict['PDB'])

    from .CG import buildENMPreset
    if params_dict['CG']['preset']=='bbENM':
        params_dict['CG']['chain_starts'] = chain_starts
    kirch, hessian = buildENMPreset(coords,  **params_dict['CG'])

    from .NMA import modeCalc
    evals, evecs = modeCalc(hessian, **params_dict['NMA'])

    from .NMA import fitCompareBfactors
    evals_scaled, evecs_scaled, cc, gamma, n_modes = fitCompareBfactors(evals, evecs, bfactors, pdb, **params_dict['b_factors'])

    from .NMA import calcDistFlucts
    from .QRC import findQuasiRigidClusters

    dist_flucts = calcDistFlucts(evals_scaled, evecs, coords)


    if not 'cluster_stop' in params_dict['QRC'].keys():
        print('No cluster range specified, defaulting to: 4-[number of chains]')
        import biotite.structure as struc
        params_dict['QRC']['cluster_start'] = 4
        params_dict['QRC']['cluster_stop'] = struc.get_chain_count(capsid)+1
        params_dict['QRC']['cluster_step'] = 2

    if 'clust_ignore_chains' in params_dict['QRC'].keys():
        from .QRC import filterChains
        clust_mask = filterChains(calphas, clust_ignore_chains=params_dict['QRC']['clust_ignore_chains'])
        params_dict['QRC'].pop('clust_ignore_chains')
    else:
        clust_mask = None

    labels, score, residue_scores = findQuasiRigidClusters(pdb, dist_flucts, clust_mask=clust_mask, **params_dict['QRC'])

    pycap_time = timer() - pycap_start
    print(f'pyCapsid total execution time for {pdb}: {pycap_time}')

    vis_method = params_dict['VIS']['method']

    from .VIS import visualizeResults

    if vis_method == 'chimerax':
        visualizeResults(pdb, capsid, labels, clust_mask=clust_mask, method=vis_method, chimerax_path=params_dict['VIS']['chimerax_path'])
    elif vis_method == 'nglview':
        from pyCapsid.VIS import createCapsidView
        view_clusters = createCapsidView(pdb, capsid)
        literal = """To visualize results in notebook use the following code in another cell:
        from pyCapsid.VIS import createClusterRepresentation
        createClusterRepresentation(pdb, labels, view_clusters)
        view_clusters
        """
        print(literal)
        return pdb, labels, view_clusters
    else:
        print('method must be chimerax or nglview')