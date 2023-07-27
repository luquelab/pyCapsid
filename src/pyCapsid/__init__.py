

def read_config(file_path):
    import toml
    params_dict = toml.load(file_path)

    return params_dict


def create_directories(params_dict):
    from pathlib import Path

    # save_params = ['save_pdb_path', 'save_cg_path', 'save_mode_path', 'save_bfactors_path', 'save_results_path']
    if 'save_pdb_path' not in params_dict['PDB'].keys():
        params_dict['PDB']['save_pdb_path'] = '/'
    if 'save_cg_path' not in params_dict['CG'].keys():
        params_dict['CG']['save_cg_path'] = '/'
    if 'save_mode_path' not in params_dict['NMA'].keys():
        params_dict['NMA']['save_mode_path'] = '/'
    if 'save_bfactors_path' not in params_dict['b_factors'].keys():
        params_dict['b_factors']['save_bfactors_path'] = '/'
    if 'save_results_path' not in params_dict['QRC'].keys():
        params_dict['QRC']['save_results_path'] = '/'

    save_suffix = '_path'
    params = params_dict.copy()

    top_path = params['PDB']['save_all_path']
    Path(top_path).mkdir(parents=True, exist_ok=True)

    # This needs refactoring but it does its job of creating all the folders for saving results

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


def run_capsid_report(params_path):
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


def chimeraxReportScriptVis(report_path, remote=True, chimerax_path=None, **kwargs):
    """Launches ChimeraX and runs a script that visualizes the results.

    :param labels:
    :param pdb:
    :param remote:
    :param chimerax_path:
    :param pdb_path:
    :param save_path:
    :param script_path:
    :param labels_path:
    :return:
    """

    import os
    import platform
    if chimerax_path is None:
        print('No chimerax path specified, checking in default locations for your OS')
        if platform.system() == 'Linux':
            chimerax_path = 'usr/bin/chimerax'
        elif platform.system() == 'Windows':
            chimerax_path = 'C:\\Program Files\\ChimeraX\\bin\\ChimeraX.exe'
        elif platform.system() == 'Darwin':
            chimerax_path = '/Applications/ChimeraX'
        else:
            print('No chimerax path is given and cannot check default locations')
            chimerax_path = ''

    # get path to chimerax script
    import os

    script_path = f'{report_path}/chimerax/chimerax_script_colab.py'

    chimerax_exe = chimerax_path  # + '\\ChimeraX.exe'
    if platform.system() == 'Windows':
        cmd_string = f'""{chimerax_exe}" --script "{script_path}""'
    elif platform.system() == 'Linux':
        cmd_string = f'{chimerax_exe} --script "{script_path}"'
    else:
        cmd_string = f'{chimerax_exe} --script "{script_path}"'
    print(cmd_string)
    os.system(cmd_string)


def run_capsid_report(params_path):
    from timeit import default_timer as timer
    pycap_start = timer()

    params_dict = read_config(params_path)
    create_directories(params_dict)
    save_all_path = params_dict['PDB'].pop('save_all_path')

    if 'suppress_plots' in params_dict['plotting'].keys():
        if params_dict['plotting']['suppress_plots']:
            import matplotlib
            matplotlib.use('Agg')

    from pyCapsid.PDB import getCapsid
    pdb = params_dict['PDB']['pdb']
    capsid, calphas, asymmetric_unit, coords, bfactors, chain_starts, title = getCapsid(**params_dict['PDB'])

    from pyCapsid.CG import buildENMPreset
    if params_dict['CG']['preset'] == 'bbENM':
        params_dict['CG']['chain_starts'] = chain_starts
    kirch, hessian = buildENMPreset(coords, **params_dict['CG'])

    from pyCapsid.NMA import modeCalc
    evals, evecs = modeCalc(hessian, **params_dict['NMA'])

    from pyCapsid.NMA import fitCompareBfactors
    evals_scaled, evecs_scaled, cc, gamma, n_modes = fitCompareBfactors(evals, evecs, bfactors, pdb,
                                                                        **params_dict['b_factors'])

    from pyCapsid.NMA import calcDistFlucts
    from pyCapsid.QRC import findQuasiRigidClusters

    dist_flucts = calcDistFlucts(evals_scaled, evecs, coords)

    if not 'cluster_stop' in params_dict['QRC'].keys():
        print('No cluster range specified, defaulting to: 4-[number of chains]')
        import biotite.structure as struc
        params_dict['QRC']['cluster_start'] = 4
        params_dict['QRC']['cluster_stop'] = struc.get_chain_count(capsid) + 1
        params_dict['QRC']['cluster_step'] = 2

    if 'clust_ignore_chains' in params_dict['QRC'].keys():
        from pyCapsid.QRC import filterChains
        clust_mask = filterChains(calphas, clust_ignore_chains=params_dict['QRC']['clust_ignore_chains'])
        params_dict['QRC'].pop('clust_ignore_chains')
    else:
        clust_mask = None

    labels, score, residue_scores = findQuasiRigidClusters(pdb, dist_flucts, clust_mask=clust_mask,
                                                           **params_dict['QRC'])

    pycap_time = timer() - pycap_start
    print(f'pyCapsid total execution time for {pdb}: {pycap_time}')

    vis_method = params_dict['VIS']['method']

    from pyCapsid.VIS import visualizeResults

    createReport(pdb, save_all_path, n_modes, residue_scores, asymmetric_unit, calphas, params_dict['CG']['preset'], cc,
                 gamma, n_modes)

    if vis_method == 'chimerax':
        report_path = f'{save_all_path}/{pdb}_pyCapsid_report'
        chimeraxReportScriptVis(report_path, chimerax_path=params_dict['VIS']['chimerax_path'])
        import webbrowser
        import os
        filename = 'file:///' + os.path.abspath(f'{report_path}/pyCapsid_report.html')
        webbrowser.open_new_tab(filename)
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


def createReport(pdb, save_all_path, n_modes_best, residue_scores, asymmetric_unit, calphas, ENM_model, cc, gamma,
                 n_modes):
    # Prepare file names
    file_name = 'pyCapsid_report'
    file_md = file_name + '.md'  # Markdown
    file_html = file_name + '.html'  # HTML
    file_docx = file_name + '.docx'

    # Include libraries to generate plots
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np

    import os
    parent_dir = save_all_path
    report_dir_name = f'{pdb}_pyCapsid_report'
    report_dir = parent_dir + '/' + report_dir_name

    if not os.path.exists(report_dir):
        os.mkdir(report_dir)

    # Generate figures folder
    figures_dir_name = 'figures'
    figures_dir = report_dir + '/' + figures_dir_name

    if not os.path.exists(figures_dir):
        os.mkdir(figures_dir)

    # Output files (quantiative plots)

    # Correlation coefficients (CC) across modes
    cc_modes_name_out = 'cc_modes'
    cc_modes_fig_png = cc_modes_name_out + '.png'
    cc_modes_fig_svg = cc_modes_name_out + '.svg'
    cc_modes_fig_eps = cc_modes_name_out + '.eps'
    cc_modes_fig_pdf = cc_modes_name_out + '.pdf'
    cc_modes_data_csv = cc_modes_name_out + '.csv'
    cc_modes_dir = figures_dir + '/' + cc_modes_name_out
    if not os.path.exists(cc_modes_dir):
        os.mkdir(cc_modes_dir)

    # b-factors
    b_factors_name_out = 'b_factors'
    b_factors_fig_png = b_factors_name_out + '.png'
    b_factors_fig_svg = b_factors_name_out + '.svg'
    b_factors_fig_eps = b_factors_name_out + '.eps'
    b_factors_fig_pdf = b_factors_name_out + '.pdf'

    b_factors_fit_name_out = 'b_factors_fit'
    b_factors_fit_fig_png = b_factors_fit_name_out + '.png'
    b_factors_fit_fig_svg = b_factors_fit_name_out + '.svg'
    b_factors_fit_fig_eps = b_factors_fit_name_out + '.eps'
    b_factors_fit_fig_pdf = b_factors_fit_name_out + '.pdf'

    b_factors_data_csv = b_factors_name_out + '.csv'
    b_factors_dir = figures_dir + '/' + b_factors_name_out
    if not os.path.exists(b_factors_dir):
        os.mkdir(b_factors_dir)

    # Cluster quality score (CQS)
    cluster_quality_name_out = 'cluster_quality'
    cluster_quality_fig_png = cluster_quality_name_out + '.png'
    cluster_quality_fig_svg = cluster_quality_name_out + '.svg'
    cluster_quality_fig_eps = cluster_quality_name_out + '.eps'
    cluster_quality_fig_pdf = cluster_quality_name_out + '.pdf'
    cluster_quality_data_csv = cluster_quality_name_out + '.csv'
    cluster_quality_dir = figures_dir + '/' + cluster_quality_name_out
    if not os.path.exists(cluster_quality_dir):
        os.mkdir(cluster_quality_dir)

    # Cluster quality score colorbar
    cluster_quality_colorbar_name_out = 'quality_score_colorbar'
    cluster_quality_colorbar_fig_svg = cluster_quality_colorbar_name_out + '.svg'
    cluster_quality_colorbar_fig_png = cluster_quality_colorbar_name_out + '.png'

    # Generate correlation coefficients dataframe and plots

    ## Upload and store file
    dir_loc = parent_dir
    file_name = 'CC_by_mode.npz'
    path_target = dir_loc + '/' + file_name
    cc_modes = np.load(path_target)

    ## Generate and save data frame
    dir_loc = cc_modes_dir
    file_name = cc_modes_data_csv
    path_target = dir_loc + '/' + file_name

    npz = cc_modes
    df = pd.DataFrame(data=[npz['mode_indices'], npz['data']]).T
    df.columns = ['Modes', 'CC']
    df_cc_modes = df
    df_cc_modes.to_csv(path_target, index=False)
    cc_best = df_cc_modes[df_cc_modes['Modes'] == n_modes_best]['CC'].values[0]

    ## Generate plots
    plt.ioff()  # Comment this command to display figure in the notebook
    sns.set_style('darkgrid')
    fig, ax = plt.subplots(1, 1)
    sns.lineplot(ax=ax, x='Modes', y='CC', data=df_cc_modes)
    ax.set_title('Agreement between predicted and empirical B-factors')
    ax.set_xlabel('Number of modes')
    ax.set_ylabel('Correlation coefficient (CC)')
    ax.vlines(x=n_modes_best, ymin=0, ymax=cc_best, color='black')
    # plt.show(fig) ## Comment out to display figure in notebook

    ## Save plots
    dir_loc = cc_modes_dir
    ### PNG
    file_name = cc_modes_fig_png
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='png', dpi=300)
    ### EPS
    file_name = cc_modes_fig_eps
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='eps', bbox_inches='tight')
    ### SVG
    file_name = cc_modes_fig_svg
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='svg')
    ### PDF
    file_name = cc_modes_fig_pdf
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='pdf')

    ## Generate caption
    text = '**Figure: Selection of the best number of modes** (below). Correlation coefficient obtained comparing the experimental B-factors with respect the B-factors predicted an increasing number of low frequency modes. The black line higlights the number of modes that yield the optimal correlation coefficient with the empirical B-factors.'
    caption_cc_modes = text

    # Generate b-factors dataframe and plots

    ## Upload and store file
    dir_loc = parent_dir
    file_name = 'b_factors.npz'
    path_target = dir_loc + '/' + file_name
    b_factors = np.load(path_target)

    ## Generate and save data frame
    ## (!! revise method for residues in asymmetric unit)
    dir_loc = b_factors_dir
    file_name = b_factors_data_csv
    path_target = dir_loc + '/' + file_name

    npz = b_factors
    k = npz['k']
    unscaled = npz['bfactors_predicted'] / k
    x_range = np.arange(unscaled.min(), unscaled.max(), 0.01)
    y_line = x_range * k
    df_bfactors = pd.DataFrame(
        data=[npz['residue_numbers'], npz['bfactors_experimental'], npz['bfactors_predicted'], unscaled, x_range,
              y_line]).T
    df_bfactors.columns = ['Residue', 'Experimental', 'Predicted', 'Unscaled', 'x_line', 'y_line']
    n_res_asym = int(len(df_bfactors) / 60)  # ! residues in asymmetric unit (refactor)
    df_bfactors_asym = df_bfactors[0:n_res_asym]
    # len(df_bfactors_asym)
    df_bfactors_asym.to_csv(path_target, index=False)

    ## Generate plots
    df = df_bfactors_asym
    plt.ioff()  # Comment this command to display figure in the notebook
    sns.set_style('darkgrid')
    fig, ax = plt.subplots(1, 1, figsize=(6, 4), sharex=False)

    # b-factor profile
    sns.lineplot(ax=ax, x='Residue', y='Predicted', data=df, label='Predicted')
    sns.lineplot(ax=ax, x='Residue', y='Experimental', data=df, label='Experimental')
    ax.set_title('Protein complex\'s vibration profile (asymmetric unit)')
    ax.set_xlabel('Residue number')
    ax.set_ylabel(r"B-factor ($\AA^2$)")
    ax.legend(frameon=False)

    # plt.show(fig) ## Comment out to display figure in notebook

    ## Save plots
    dir_loc = b_factors_dir
    ### PNG
    file_name = b_factors_fig_png
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='png', dpi=300)
    ### EPS
    file_name = b_factors_fig_eps
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='eps', bbox_inches='tight')
    ### SVG
    file_name = b_factors_fig_svg
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='svg')
    ### PDF
    file_name = b_factors_fig_pdf
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='pdf')

    ## Generate caption
    text = '**Figure: Correlation of B-factors** (below). Empirical (blue) and predicted (orange) B-factors for each residue in the asymmetric unit.'
    caption_bfactors = text

    # Generate b-factors dataframe and plots

    ## Upload and store file
    dir_loc = parent_dir
    file_name = 'b_factors.npz'
    path_target = dir_loc + '/' + file_name
    b_factors = np.load(path_target)

    ## Generate and save data frame
    ## (!! revise method for residues in asymmetric unit)
    dir_loc = b_factors_dir
    file_name = b_factors_data_csv
    path_target = dir_loc + '/' + file_name

    npz = b_factors
    k = npz['k']
    unscaled = npz['bfactors_predicted'] / k
    x_range = np.arange(unscaled.min(), unscaled.max(), 0.01)
    y_line = x_range * k
    df_bfactors = pd.DataFrame(
        data=[npz['residue_numbers'], npz['bfactors_experimental'], npz['bfactors_predicted'], unscaled, x_range,
              y_line]).T
    df_bfactors.columns = ['Residue', 'Experimental', 'Predicted', 'Unscaled', 'x_line', 'y_line']
    n_res_asym = int(len(df_bfactors) / 60)  # ! residues in asymmetric unit (refactor)
    df_bfactors_asym = df_bfactors[0:n_res_asym]
    # len(df_bfactors_asym)
    df_bfactors_asym.to_csv(path_target, index=False)

    ## Generate plots
    df = df_bfactors_asym
    plt.ioff()  # Comment this command to display figure in the notebook
    sns.set_style('darkgrid')
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

    # b-factor fit
    sns.scatterplot(ax=ax, x='Unscaled', y='Experimental', data=df, label='Scatter')
    sns.lineplot(ax=ax, x='x_line', y='y_line', data=df, label='Linear Fit')
    ax.set_title('ENM Calibration (linear fit)')
    ax.set_xlabel('Predicted squared fluctuations for $\gamma  = 1$ ($\AA^2$)')
    ax.set_ylabel(r"B-factor ($\AA^2$)")
    ax.legend(frameon=False)

    # plt.show(fig) ## Comment out to display figure in notebook

    ## Save plots
    dir_loc = b_factors_dir
    ### PNG
    file_name = b_factors_fit_fig_png
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='png', dpi=300)
    ### EPS
    file_name = b_factors_fit_fig_eps
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='eps', bbox_inches='tight')
    ### SVG
    file_name = b_factors_fit_fig_svg
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='svg')
    ### PDF
    file_name = b_factors_fit_fig_pdf
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='pdf')

    ## Generate caption
    text = '**Figure: Model calibration** (below). The empirical B-factors are plotted against the predicted fluctuations. A regression line was fitted. From the slope (a), the calibration constant (gamma) was obtained using the formula: *gamma = (8(pi)^2)/(3a)*.'
    caption_bfactors_fit = text

    # Generate quality score dataframe and plots

    ## Obtain optimal number of clusters
    dir_loc = parent_dir
    file_name = str(pdb) + '_final_results.npz'
    path_target = dir_loc + '/' + file_name
    npz_loaded = np.load(path_target)
    nc = int(npz_loaded['nc'])

    ## Upload and store quality score file
    dir_loc = parent_dir
    file_name = str(pdb) + '_final_results_full.npz'
    path_target = dir_loc + '/' + file_name
    npz_loaded = np.load(path_target)

    ## Generate and save data frame
    dir_loc = cluster_quality_dir
    file_name = cluster_quality_data_csv
    path_target = dir_loc + '/' + file_name

    df = pd.DataFrame(data=[npz_loaded['nc_range'], npz_loaded['score'], npz_loaded['numtypes']]).T
    df.columns = ['Clusters', 'Quality score', 'Unique']
    df.to_csv(path_target, index=False)

    # Obtain score and unique clustering for optimal number of clusters (nc)
    ii = df[df['Clusters'] == nc].index[0]
    score = df['Quality score'].iloc[ii]
    unique = df['Unique'].iloc[ii]

    ## Generate and save plots
    plt.ioff()  # Comment this command to display figure in the notebook
    sns.set_style('darkgrid')
    fig, ax = plt.subplots(2, 1, figsize=(6, 3), sharex=True)
    fig.suptitle('Cluster quality score analysis')
    sns.lineplot(ax=ax[0], x='Clusters', y='Quality score', data=df)
    sns.lineplot(ax=ax[1], x='Clusters', y='Unique', data=df)
    ax[1].set_xlabel('Number of clusters')
    ax[0].set_ylabel('Quality score')
    ax[1].set_ylabel('Unique clusters')
    ax[0].vlines(x=nc, ymin=0, ymax=score, color='black')
    ymax = ax[1].get_ylim()[1]
    ax[1].vlines(x=nc, ymin=unique, ymax=ymax, color='black')
    fig.tight_layout()
    # plt.show(fig) ## Comment out to display figure in notebook

    ## Save plots
    dir_loc = cluster_quality_dir
    ### PNG
    file_name = cluster_quality_fig_png
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='png', dpi=300)
    ### EPS
    file_name = cluster_quality_fig_eps
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='eps', bbox_inches='tight')
    ### SVG
    file_name = cluster_quality_fig_svg
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='svg')
    ### PDF
    file_name = cluster_quality_fig_pdf
    path_target = dir_loc + '/' + file_name
    fig.savefig(path_target, format='pdf')

    ## Generate caption
    text = '**Figure: Optimal clustering selection** (below). Quality score (top) and number of unique clusters (bottom) obtained as a function of the number of clusters tested. The black lines in each plot higlight the quality score and number of unique clusters for the optimal number of clusters identified.'
    caption_cluster_quality = text

    # Quality score colorbar figure
    import matplotlib.colorbar as colorbar
    import matplotlib.pyplot as plt
    from pyCapsid.VIS import clusters_colormap_hexcolor
    import numpy as np
    hexcolor, cmap = clusters_colormap_hexcolor(residue_scores, rwb_scale=True)
    fig, ax = plt.subplots(figsize=(10, 1))
    cb = colorbar.ColorbarBase(ax, orientation='horizontal',
                               cmap=cmap, norm=plt.Normalize(np.min(residue_scores), np.max(residue_scores)))
    # plt.show()
    fig.suptitle('Clustering quality scores ranging from low (red) to high (blue)', y=-0.02)
    fig.tight_layout()
    fig.savefig(cluster_quality_colorbar_fig_svg)
    fig.savefig(f'{cluster_quality_dir}/{cluster_quality_colorbar_fig_svg}')

    # Structures folder
    # Generate folder
    structures_dir = figures_dir + '/' + 'structures'
    if not os.path.exists(structures_dir):
        os.mkdir(structures_dir)

    # Set structure view
    asym_unit_name_out = f'{pdb}_asymmetric_unit'
    asym_unit_fig_png = asym_unit_name_out + '.png'

    # Save image
    dir_loc = structures_dir
    file_name = asym_unit_fig_png
    path_target = dir_loc + '/' + file_name
    from PIL import Image
    image = Image.new('RGB', (20, 20))
    image.save(path_target)

    # Set structure view
    capsid_name_out = f'{pdb}_full_capsid'
    capsid_fig_png = capsid_name_out + '.png'

    # Save image
    dir_loc = structures_dir
    file_name = capsid_fig_png
    path_target = dir_loc + '/' + file_name

    from PIL import Image
    image = Image.new('RGB', (20, 20))
    image.save(path_target)

    # Set structure view
    capsid_scores_name_out = f'{pdb}_residue_cluster_scores'
    capsid_scores_fig_png = capsid_scores_name_out + '.png'

    # Save image
    dir_loc = structures_dir
    file_name = capsid_scores_fig_png
    path_target = dir_loc + '/' + file_name

    from PIL import Image
    image = Image.new('RGB', (20, 20))
    image.save(path_target)

    capsid_clusters_name_out = f'{pdb}_highest_quality_clusters'
    capsid_clusters_fig_png = capsid_clusters_name_out + '.png'

    # Save image
    dir_loc = structures_dir
    file_name = capsid_clusters_fig_png
    path_target = dir_loc + '/' + file_name

    from PIL import Image
    image = Image.new('RGB', (20, 20))
    image.save(path_target)

    chimerax_dir = report_dir + "/chimerax/"
    if not os.path.exists(chimerax_dir):
        os.mkdir(chimerax_dir)

    import pyCapsid.scripts as scpath
    import os
    import shutil

    script_path = os.path.abspath(scpath.__file__)
    script_path_1 = script_path.replace('__init__', 'chimerax_script_colab')
    script_path_2 = script_path.replace('__init__', 'chimerax_script_animate_mode')
    shutil.copy(script_path_1, chimerax_dir)
    shutil.copy(script_path_2, chimerax_dir)
    shutil.copy(f'{save_all_path}/{pdb}_final_results_full.npz', chimerax_dir)
    shutil.copy(f'{save_all_path}/{pdb}_final_results.npz', chimerax_dir)
    # shutil.copy('modes.npz', chimerax_dir + f'{pdb}_modes.npz')
    # if pdb_source == 'upload':
    #   shutil.copy(pdb_file_name, chimerax_dir)

    # Include only the first non-degenerate eigenmode in the chimeras folder
    modes_saved = np.load(f'{save_all_path}/modes.npz')
    evals_saved = modes_saved['eigen_vals']
    evecs_saved = modes_saved['eigen_vecs']
    uniques, inds, counts = np.unique(evals_saved.round(decimals=8), return_index=True, return_counts=True)
    icoEvalInds = inds[counts == 1]
    print('Indices of non-degenerate eigenmodes: ', icoEvalInds)

    if len(icoEvalInds) == 0:
        print('No non-degenerate modes found, including lowest-frequency mode for visualization')
        n_mode = 0
    else:
        print('Including modes up to lowest frequency non-degenerate mode for visualization')
        n_mode = icoEvalInds[0] + 1

    np.savez_compressed(chimerax_dir + f'{pdb}_modes.npz', eigen_vals=evals_saved[:n_mode],
                        eigen_vecs=evecs_saved[:, :n_mode])

    # Write readme
    with open(chimerax_dir + "readme_chimerax.md", 'w') as file:
        line = "# ChimeraX visualization script"
        file.write('\n' + line)
        line = 'To visualize the clustering results in ChimeraX we provide a python script in the pyCapsid report folder. To use this script, after extracting the pyCapsid_report folder, open ChimeraX and run the following command:'
        file.write('\n' + line)
        line = f'``` \n runscript "path/to/pyCapsid_report/chimerax/chimerax_script_colab.py" \n ```'
        file.write('\n' + line)
        line = f'``` Alternatively you can use the following command, which will open a window and  prompt you to find and select the script:'
        file.write('\n' + line)
        line = f'``` \n runscript browse \n ```'
        file.write('\n' + line)
        line = 'This will create the 4 images in the pyCapsid_report/figures/structures/ folder which will be displayed in the report.'
        file.write('\n' + line)
        line = 'The user may also visualize motion of the structure along a given normal mode using a separate script in a similar manner:'
        file.write('\n' + line)
        line = f'``` \n runscript "path/to/pyCapsid_report/chimerax/chimerax_script_animate_mode.py" \n ```'
        file.write('\n' + line)
        line = f'``` \n runscript browse \n ```'
        file.write('\n' + line)
        line = f'This will create an .mp4 file in the figures/structures folder of the report.'
        file.write('\n' + line)
        line = f'Only one mode may be visualized at a time. By default the lowest frequency non-degenerate mode is visualized, and only modes of that frequency and lower are included in the report. To visualize higher frequency modes, you may copy the `modes.npz` file from `pyCapsid_results` to the `chimerax` folder of the report and rename it accordingly.'
        file.write('\n' + line)
        line = f'Both scripts can be modified to fit the users needs. To do this reference the ChimeraX commands (https://www.cgl.ucsf.edu/chimerax/docs/user/index.html) and the ChimeraX Python API (https://www.cgl.ucsf.edu/chimerax/docs/devel/index.html#).'
        file.write('\n' + line)

    # Generate markdown output file
    dir_loc = report_dir
    file_name = file_md
    path_target = dir_loc + '/' + file_name
    f = open(path_target, 'w')

    # Load data needed
    from datetime import date
    today = date.today()

    # HEADING
    text = '# pyCapsid Report'
    f.write('\n' + text)

    text = today.strftime("%B %d, %Y")
    f.write('\n' + text)

    # INPUT STRUCTURE
    f.write('\n')
    text = '## Input structure'
    f.write('\n' + text)

    text = 'Identifier: ' + str(pdb)
    f.write('\n' + text)

    f.write('\n')
    import biotite.structure as struc
    n_asym_residues = struc.get_residue_count(asymmetric_unit)
    text = 'Number of residues in the asymmetric unit: ' + str(n_asym_residues)
    f.write('\n' + text)

    f.write('\n')
    n_asym_chains = struc.get_chain_count(asymmetric_unit)
    text = 'Number of protein chains in the asymmetric unit: ' + str(n_asym_chains)
    f.write('\n' + text)

    f.write('\n')
    n_full_residues = len(calphas)
    asym_factor = int(n_full_residues / n_asym_residues)
    text = 'Multiplying factor to generate the complete protein complex: ' + str(asym_factor)
    f.write('\n' + text)

    f.write('\n')
    text = '+ If the multiplying factor is one (*m = 1*), the protein complex and the asymmetric unit are the same.'
    f.write('\n' + text)

    f.write('\n')
    text = '+ If the multiplying factor is larger than one (*m > 1*), the protein complex is *m* times the asymmetric unit.'
    f.write('\n' + text)

    f.write('\n')
    text = 'Number of residues in the protein complex: ' + str(n_full_residues)
    f.write('\n' + text)

    f.write('\n')
    n_full_chains = int(asym_factor * n_asym_chains)
    text = 'Number of protein chains in the protein complex: ' + str(n_full_chains)
    f.write('\n' + text)

    ## Asymmetric unit
    f.write('\n')
    text = '***'
    f.write('\n' + text)
    f.write('\n')
    text = '**Figure: Asymmetric unit** (below). Ribbon diagram of the protein complex\'s asymmetric unit.'
    f.write('\n' + text)
    f.write('\n')
    text = '![Asymmetric unit](' + './figures/structures/' + asym_unit_fig_png + ')'
    f.write('\n' + text)

    ## Full structure
    f.write('\n')
    text = '***'
    f.write('\n' + text)
    f.write('\n')
    text = '**Figure: Full protein complex** (below). Ribbon diagram of the full protein complex.'
    f.write('\n' + text)
    f.write('\n')
    text = '![Full protein complex](' + './figures/structures/' + capsid_fig_png + ')'
    f.write('\n' + text)
    f.write('\n')
    text = '***'
    f.write('\n' + text)

    # Elastic network model
    f.write('\n')
    text = '## Elastic network model'
    f.write('\n' + text)

    f.write('\n')
    text = 'Elastic model used: ' + str(ENM_model)
    f.write('\n' + text)

    f.write('\n')
    text = 'Calibrated stiffness constant (gamma): ' + str(round(gamma, 2))
    f.write('\n' + text)
    f.write('\n')
    text = '+ This constant was fitted to scale the model to the structure, assuming a linear relationship between the residues fluctuations and B-factors.'
    f.write('\n' + text)

    f.write('\n')
    text = '***'
    f.write('\n' + text)
    f.write('\n')
    text = caption_bfactors_fit
    f.write('\n' + text)
    f.write('\n')
    ### Path to render figure by user locally in markdown
    subdir_sep = b_factors_dir.split(sep='/')[::-1][0:2][::-1]
    subdir = './' + subdir_sep[0] + '/' + subdir_sep[1]
    text = '![Model calibration](' + subdir + '/' + b_factors_fit_fig_svg + ')'
    f.write('\n' + text)
    ### Path to render figure by pypandoc in Colab
    dir_loc = b_factors_dir
    file_name = b_factors_fit_fig_svg
    path_target = dir_loc + '/' + file_name
    text = '![](' + path_target + ')'
    f.write('\n' + text)

    # Normal modes and B-factors
    CC = float(b_factors['CC'])

    f.write('\n')
    text = '## Normal mode analysis (NMA)'
    f.write('\n' + text)

    f.write('\n')
    text = 'Optimal number of modes reproducing B-factors: ' + str(n_modes_best)
    f.write('\n' + text)

    f.write('\n')
    text = 'Correlation between empirical and predicted B-factors: ' + str(round(CC, 2))
    f.write('\n' + text)

    ## CC vs number of modes figure
    f.write('\n')
    text = '***'
    f.write('\n' + text)
    f.write('\n')
    text = caption_cc_modes
    f.write('\n' + text)
    f.write('\n')
    ### Path to render figure by user locally in markdown
    subdir_sep = cc_modes_dir.split(sep='/')[::-1][0:2][::-1]
    subdir = './' + subdir_sep[0] + '/' + subdir_sep[1]
    text = '![Selection of the best number of modes](' + subdir + '/' + cc_modes_fig_svg + ')'
    f.write('\n' + text)
    ### Path to render figure by pypandoc in Colab
    dir_loc = cc_modes_dir
    file_name = cc_modes_fig_svg
    path_target = dir_loc + '/' + file_name
    text = '![](' + path_target + ')'
    f.write('\n' + text)

    ## B-factors figure
    f.write('\n')
    text = caption_bfactors
    f.write('\n' + text)

    f.write('\n')
    ### Path to render figure by user locally in markdown
    subdir_sep = b_factors_dir.split(sep='/')[::-1][0:2][::-1]
    subdir = './' + subdir_sep[0] + '/' + subdir_sep[1]
    text = '![Correlation of B-factors](' + subdir + '/' + b_factors_fig_svg + ')'
    f.write('\n' + text)
    ### Path to render figure by pypandoc in Colab
    dir_loc = b_factors_dir
    file_name = b_factors_fig_svg
    path_target = dir_loc + '/' + file_name
    text = '![](' + path_target + ')'
    f.write('\n' + text)

    f.write('\n')
    text = '***'
    f.write('\n' + text)

    # Quasi-rigid units
    ## Upload and store file
    dir_loc = parent_dir
    file_name = str(pdb) + '_final_results.npz'
    path_target = dir_loc + '/' + file_name
    npz_loaded = np.load(path_target)
    nc = int(npz_loaded['nc'])

    f.write('\n')
    text = '## Quasi-rigid mechanical units'
    f.write('\n' + text)

    f.write('\n')
    text = 'Number of optimal quasi-rigid mechanical units identified: ' + str(nc)
    f.write('\n' + text)

    # f.write('\n')
    # text = '!!! Add figure plotting the quality score and number of unique clusters as a function of the number of selected clusters. Include captions following standard publication style.'
    # f.write('\n'+text)

    ## Cluster quality figure
    f.write('\n')

    f.write('\n')
    text = caption_cluster_quality
    f.write('\n' + text)

    ### Path to render figure by user locally in markdown
    subdir_sep = cluster_quality_dir.split(sep='/')[::-1][0:2][::-1]
    subdir = './' + subdir_sep[0] + '/' + subdir_sep[1]
    text = '![](' + subdir + '/' + cluster_quality_fig_svg + ')'
    f.write('\n' + text)
    ### Path to render figure by pypandoc in Colab
    dir_loc = cluster_quality_dir
    file_name = cluster_quality_fig_svg
    path_target = dir_loc + '/' + file_name
    text = '![](' + path_target + ')'
    f.write('\n' + text)

    f.write('\n')
    text = '***'
    f.write('\n' + text)
    f.write('\n')
    text = '**Figure: Quasi-rigid clusters** (below). Ribbon representation of the complete structure, with each residue colored according to its cluster membership. Residues with the same color are members of the same quasi-rigid cluster.'
    f.write('\n' + text)
    f.write('\n')
    text = '![](' + './figures/structures/' + capsid_clusters_fig_png + ')'
    f.write('\n' + text)

    # f.write('\n')
    # text = '!!! Add figure displaying a representative of each unique cluster. Include caption following standard publication style.'
    # f.write('\n'+text)

    f.write('\n')
    text = '**Figure: Residue cluster quality score** (below). Ribbon representation of the complete structure, with each residue colored according to its cluster quality score.  This is a measure of how rigid each residue is with respect to its cluster. Blue residues make up the cores of rigid clusters, and red residues represent borders between clusters.'
    f.write('\n' + text)

    ## clusters
    f.write('\n')
    text = '![](' + './figures/structures/' + capsid_scores_fig_png + ')'
    f.write('\n' + text)

    # colorbar
    f.write('\n')
    text = '![](' + './figures/cluster_quality/' + cluster_quality_colorbar_fig_svg + ')'
    f.write('\n' + text)

    # f.write('\n')
    # text = '!!! Add figure displaying the quality score of each residue in the representative unique clusters. Include caption following standard publication style.'
    # f.write('\n'+text)

    f.write('\n')
    text = "# Missing images in the report? \n "
    f.write('\n' + text)

    f.write('\n')
    text = "To render the missing images locally using ChimeraX, check the 'chimerax' folder of the report, and follow the instructions provided in the readme file. This will automatically place the images in the .md and .html versions. For the word document, copy and paste the relevant images directly into the document."
    f.write('\n' + text)

    f.write('\n')
    text = "If you were visualizing the results using NGLView, you should have several images downloaded corresponding to the relevant images. After extracting \"pyCapsid_report.zip\", place the corresponding images from your download folder into \"pyCapsid_report/figures/structures/\" replacing the empty images in that folder. This will fix the issue for the .md and .html versions. For the word document, copy and paste the relevant images directly into the document."
    f.write('\n' + text)

    # Close ouptut file
    f.close()

    import markdown
    # Generate HTML version
    dir_loc = report_dir
    file_name = file_md
    path_in = dir_loc + '/' + file_name
    file_name = file_html
    path_out = dir_loc + '/' + file_name

    markdown.markdownFromFile(input=path_in, output=path_out)