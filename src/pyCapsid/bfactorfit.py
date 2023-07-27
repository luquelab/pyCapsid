import numpy as np
import numba as nb


# @nb.njit()
def fastFlucts(evals, evecs, n_modes, is3d):
    """

    :param evals:
    :param evecs:
    :param n_modes:
    :param is3d:
    :return:
    """
    n_d = evecs.shape[0]
    flucts = np.zeros(n_d)
    for i in range(n_modes):
        flucts += 1 / evals[i] * evecs[:, i] ** 2
    if is3d:
        return np.reshape(flucts, (-1, 3)).sum(axis=-1)
    else:
        return flucts


# @nb.njit()
def checkIcoFlucts(flucts):
    """

    :param flucts:
    :return:
    """
    F = np.reshape(flucts, (60, -1))
    devs = np.ptp(F, axis=0)
    d = np.max(devs)
    # print('Maximum deviation from icosahedral: ', np.max(devs))
    return d


def springFit(bfactors, sqFlucts):
    """

    :param bfactors:
    :param sqFlucts:
    :return:
    """
    import statsmodels.api as sm

    intercept = False
    if intercept:
        sqFlucts = sm.add_constant(sqFlucts)
    M = sm.RLM(bfactors, sqFlucts, M=sm.robust.norms.HuberT())
    results = M.fit()
    # print(results.summary(xname=['Squared Fluctuations'], yname='B-factors', title='Results of fitting predicted squared fluctuations'))
    print(results.summary2(xname=['Squared Fluctuations'], yname='B-factors',
                          title='Summary of regression results: predicted squared fluctuations vs B-factors'))
    a = results.params[-1]
    stderr = results.bse * np.sqrt(bfactors.shape[0])
    pv = results.pvalues
    ci = results.conf_int(alpha=0.05)

    if intercept:
        b = results.params[0]
    else:
        b = 0

    # from sklearn.linear_model import HuberRegressor
    #
    # huber = HuberRegressor(fit_intercept=False, tol=0, alpha=0.0).fit(sqFlucts, bfactors)
    # a = huber.coef_
    # b = huber.intercept_

    return a, b, stderr, ci, pv


def fluctModes(evals, evecs, bfactors, is3d, isIco):
    """

    :param evals:
    :param evecs:
    :param bfactors:
    :param is3d:
    :return:
    """
    coeffs = []
    ico_devs = []

    for n_modes in range(1, len(evals)):
        flucts = fastFlucts(evals, evecs, n_modes, is3d)
        cc = np.corrcoef(bfactors, flucts)[1, 0]
        if isIco:
            icodev = checkIcoFlucts(flucts)
        else:
            icodev = 0
        coeffs.append(cc)
        ico_devs.append(icodev)
    return coeffs, ico_devs


def fitBfactors(evals, evecs, bfactors, is3d, isIco=True, fitModes=False, plotModes=False, forceIco=False, icotol=0.002, save_path='./'):
    """

    :param evals:
    :param evecs:
    :param bfactors:
    :param is3d:
    :param fitModes:
    :param plotModes:
    :param forceIco:
    :param icotol:
    :return:
    """
    n_modes = evals.shape[0]
    if plotModes:
        coeffs, ico_devs = fluctModes(evals, evecs, bfactors, is3d, isIco)
        plotByMode(np.arange(1, n_modes), coeffs, 'CC', save_path=save_path)
        if fitModes:
            if forceIco:
                plotByMode(np.arange(1, n_modes), ico_devs, 'Icosahedral Deviation', save_path=save_path)
                icoI = np.nonzero(np.array(ico_devs) < icotol)
                n_m = np.argmax(np.array(coeffs)[icoI])
                coeff = np.array(coeffs)[icoI][n_m]
                ico_dev = np.array(ico_devs)[icoI][n_m]
            else:
                n_m = np.argmax(coeffs)
                coeff = coeffs[n_m]
                ico_dev = ico_devs[n_m]
        else:
            n_m = n_modes
            flucts = fastFlucts(evals, evecs, n_m, is3d)
            coeff = np.corrcoef(bfactors, flucts)[1, 0]
            if isIco:
                ico_dev = checkIcoFlucts(flucts)
            else:
                ico_dev = 0
        flucts = fastFlucts(evals, evecs, n_m, is3d)
    else:
        n_m = n_modes
        flucts = fastFlucts(evals, evecs, n_m, is3d)
        coeff = np.corrcoef(bfactors, flucts)[1, 0]
        if isIco:
            ico_dev = checkIcoFlucts(flucts)
        else:
            ico_dev = 0

    k, intercept, stderr, ci, pv = springFit(bfactors, flucts[:, np.newaxis])

    bfactors_predicted = k * flucts + intercept

    return coeff, k, intercept, bfactors_predicted, ci, pv, ico_dev, n_m


def plotByMode(mode_indices, data, datalabel, save_path='./'):
    """

    :param mode_indices:
    :param data:
    :param datalabel:
    """
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(1, 1)
    ax.plot(mode_indices, data)
    # ax[0].vlines(nModes, np.min(coeffs), np.max(coeffs))
    ax.set_xlabel('Number Of Modes')
    ax.set_ylabel(datalabel)
    fig.suptitle(datalabel + ' vs number of low frequency modes')
    plt.show()
    fig.savefig(f'{save_path}{datalabel}_by_mode.svg')
    np.savez(f'{save_path}{datalabel}_by_mode.npz', mode_indices = mode_indices, data = data)


def fitPlotBfactors(evals, evecs, bfactors, pdb, is3d=True, fitModes=True, plotModes=False, forceIco=True, icotol=0.002,
                    save=True, save_path='./', isIco=False):
    """

    :param isIco:
    :param evals:
    :param evecs:
    :param bfactors:
    :param pdb:
    :param is3d:
    :param fitModes:
    :param plotModes:
    :param forceIco:
    :param icotol:
    :return:
    """
    coeff, k, intercept, bfactors_predicted, ci, pv, ico_dev, nmodes = fitBfactors(evals, evecs, bfactors, is3d, isIco,
                                                                                   fitModes, plotModes, forceIco,
                                                                                   icotol, save_path=save_path)

    print(f'Number of low-frequency modes used: {nmodes}')
    k_ci = np.abs(ci[0][0] - ci[0][1])
    print(f'Scale factor k between predicted fluctuations and B-factors: {k:.2e}±{k_ci:.2e}')
    gamma = (8 * np.pi ** 2) / k
    ci = (8 * np.pi ** 2) / ci

    if is3d:
        gamma = gamma / 3
        ci = ci / 3

    gamma_ci = np.abs(ci[0][0] - ci[0][1])
    print(f'Estimated spring constant gamma of ENM springs: {gamma:.2e}±{gamma_ci:.2e}')

    file = f'{save_path}b_factors.npz'
    print('Saving B-factor results in' + file)
    residue_numbers = np.arange(1, bfactors.shape[0] + 1)
    np.savez_compressed(file, bfactors_predicted=bfactors_predicted, bfactors_experimental=bfactors, CC=coeff, gamma=gamma, gamma_ci=gamma_ci, k=k, k_ci=k_ci, n_modes=nmodes, n_asym=int(bfactors.shape[0] / 60), residue_numbers=residue_numbers)
    
    import matplotlib.pyplot as plt
    import matplotlib
    fig, ax = plt.subplots(1, 1, figsize=(6, 3))
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 9}
    matplotlib.rc('font', **font)

    if isIco:
        n_asym = int(bfactors.shape[0] / 60)
        plotx_experimental = np.arange(bfactors.shape[0])[:n_asym]
        ploty_experimental = bfactors[:n_asym]
        plotx_predicted = np.arange(bfactors_predicted.shape[0])[:n_asym]
        ploty_predicted = bfactors_predicted[:n_asym]
    else:
        plotx_experimental = np.arange(bfactors.shape[0])
        ploty_experimental = bfactors
        plotx_predicted = np.arange(bfactors_predicted.shape[0])
        ploty_predicted = bfactors_predicted

    ax.plot(plotx_experimental, ploty_experimental, label='b-factors (Experimental)')
    ax.plot(plotx_predicted, ploty_predicted, label='b-factors (Predicted)')
    ax.set_ylabel(r'$Å^{2}$', fontsize=9)
    ax.set_xlabel('Residue Number', fontsize=9)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)

    ax.legend()
    fig.suptitle(
        (fr"Experimental vs Predicted b-factors: ({pdb})" + '\n' + fr"$\gamma = $ {gamma:.2e} $\pm$ {gamma_ci:.2e} "
                                                                   fr"$\frac{{k_b T}}{{Å^{2}}}$ CC = {coeff:.2f} Icosahedral deviation = {ico_dev:.2f}"),
        fontsize=9)
    fig.tight_layout()
    if save:
        fig.savefig(f'{save_path}b_factors.svg')
    plt.show()

    return coeff, gamma, nmodes
