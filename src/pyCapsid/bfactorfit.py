import numpy as np
import numba as nb


# @nb.njit()
def fastFlucts(evals, evecs, n_modes, is3d):
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
    F = np.reshape(flucts, (60, -1))
    devs = np.ptp(F, axis=0)
    d = np.max(devs)
    # print('Maximum deviation from icosahedral: ', np.max(devs))
    return d


def springFit(bfactors, sqFlucts):
    import statsmodels.api as sm

    intercept = False
    if intercept:
        sqFlucts = sm.add_constant(sqFlucts)
    M = sm.RLM(bfactors, sqFlucts, M=sm.robust.norms.HuberT())
    results = M.fit()
    print(results.summary())
    a = results.params[-1]
    stderr = results.bse * np.sqrt(bfactors.shape[0])
    pv = results.pvalues
    ci = results.conf_int(alpha=0.1)

    if intercept:
        b = results.params[0]
    else:
        b = 0

    return a, b, stderr, ci, pv


def fluctModes(evals, evecs, bfactors, is3d):
    coeffs = []
    ico_devs = []

    for n_modes in range(1, len(evals)):
        flucts = fastFlucts(evals, evecs, n_modes, is3d)
        cc = np.corrcoef(bfactors, flucts)[1, 0]
        icodev = checkIcoFlucts(flucts)
        coeffs.append(cc)
        ico_devs.append(icodev)
    return coeffs, ico_devs


def fitBfactors(evals, evecs, bfactors, is3d, fitModes=False, plotModes=False, forceIco=False, icotol=0.002):
    n_modes = evals.shape[0]
    if fitModes:
        coeffs, ico_devs = fluctModes(evals, evecs, bfactors, is3d)
        if plotModes:
            plotByMode(np.arange(1, n_modes), coeffs, 'CC')
            plotByMode(np.arange(1, n_modes), ico_devs, 'Icosahedral Deviation')

        if forceIco:
            icoI = np.nonzero(np.array(ico_devs) < icotol)
            n_m = np.argmax(np.array(coeffs)[icoI])
            coeff = np.array(coeffs)[icoI][n_m]
            ico_dev = np.array(ico_devs)[icoI][n_m]
        else:
            n_m = np.argmax(coeffs)
            coeff = coeffs[n_m]
            ico_dev = ico_devs[n_m]
        flucts = fastFlucts(evals, evecs, n_m, is3d)
    else:
        n_m = n_modes
        flucts = fastFlucts(evals, evecs, n_m, is3d)
        coeff = np.corrcoef(bfactors, flucts)[1, 0]
        ico_dev = checkIcoFlucts(flucts)

    k, intercept, stderr, ci, pv = springFit(bfactors, flucts[:, np.newaxis])

    bfactors_predicted = k * flucts + intercept

    return coeff, k, intercept, bfactors_predicted, ci, pv, ico_dev, n_m


def plotByMode(mode_indices, data, datalabel):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    ax.plot(mode_indices, data)
    # ax[0].vlines(nModes, np.min(coeffs), np.max(coeffs))
    ax.set_xlabel('Number Of Modes')
    ax.set_ylabel(datalabel)
    fig.suptitle(datalabel + ' vs number of low frequency modes')
    plt.show()


def plotBfactors(evals, evecs, bfactors, pdb, is3d=True, fitModes=False, plotModes=False, forceIco=False, icotol=0.002):
    coeff, k, intercept, bfactors_predicted, ci, pv, ico_dev, nmodes = fitBfactors(evals, evecs, bfactors, is3d,
                                                                                   fitModes,
                                                                                   plotModes, forceIco, icotol)
    print(ci)
    ci = np.abs(ci[0][0] - ci[0][1])
    gamma = (8 * np.pi ** 2) / k

    if is3d:
        gamma = gamma / 3

    import matplotlib.pyplot as plt
    import matplotlib
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 11}
    matplotlib.rc('font', **font)

    n_asym = int(bfactors.shape[0] / 60)

    ax.plot(np.arange(bfactors.shape[0])[:n_asym], bfactors[:n_asym], label='b-factors (Experimental)')
    ax.plot(np.arange(bfactors_predicted.shape[0])[:n_asym], bfactors_predicted[:n_asym], label='b-factors (Predicted)')
    ax.set_ylabel(r'$Å^{2}$', fontsize=12)
    ax.set_xlabel('Residue Number', fontsize=12)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)

    ax.legend()
    fig.suptitle(
        (fr"Experimental vs Predicted b-factors: ({pdb})" + '\n' + fr"$\gamma = $ {gamma:.2e} $\pm$ {ci:.2e} "
    fr"$\frac{{k_b T}}{{Å^{2}}}$ CC = {coeff:.2f} Icosahedral deviation = {ico_dev:.2f}"), fontsize=12)
    plt.show()

    return coeff, k, intercept, bfactors_predicted, ci, pv, ico_dev, nmodes
