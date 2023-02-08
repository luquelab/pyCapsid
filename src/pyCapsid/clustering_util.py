import numpy as np
import numba as nb


@nb.njit()
def calcCentroids(X, labels, n_clusters):
    """

    :param X:
    :param labels:
    :param n_clusters:
    :return:
    """
    centroids = np.zeros((labels.shape[0], n_clusters))
    n = n_clusters
    for i in range(n_clusters):
        mask = (labels == i)
        if not np.any(mask):
            n += -1
            centroids[i, :] = np.random.rand(n_clusters)
        else:
            clust = X[mask, :]
            cent = np.mean(clust, axis=0)
            centroids[i, :] = cent

    # if n != n_clusters:
    #     print('Some clusters unassigned')
    #     print('Assigned Clusters: ', n)

    return np.array(centroids)


@nb.njit()
def calcCosCentroids(X, labels, n_clusters):
    """

    :param X:
    :param labels:
    :param n_clusters:
    :return:
    """
    centroids = np.zeros((n_clusters, n_clusters))
    n = n_clusters
    for i in range(n_clusters):
        mask = (labels == i)
        if not np.any(mask):
            n += -1
            centroids[i, :] = np.random.rand(n_clusters)
        else:
            clust = X[mask, :]
            c = np.sum(clust, axis=0)
            cent = c / np.linalg.norm(c)
            centroids[i, :] = cent

    # if n != n_clusters:
    #     print('Some clusters unassigned')
    #     print('Assigned Clusters: ', n)

    return centroids


def median_score(coords, centroids, score_method):
    """

    :param coords:
    :param centroids:
    :param score_method:
    :return:
    """
    from sklearn.metrics import pairwise_distances
    from sklearn.preprocessing import normalize

    # n_rand = 10
    # for i in range(n_rand):
    #     randcoords = np.random.normal(size=coords.shape)
    #     randcoords = normalize(randcoords)
    #     dists = pairwise_distances(randcoords, centroids, metric='cosine')
    #     d2min = np.partition(dists, kth=2)[:, :2]
    #     if scoreMethod == 'median':
    #         b = np.median(d2min[:, 1])
    #         a = np.median(d2min[:, 0])
    #     else:
    #         b = np.mean(d2min[:, 1])
    #         a = np.mean(d2min[:, 0])
    #     r_score = b / a

    dists = pairwise_distances(coords, centroids, metric='cosine')
    cdist = pairwise_distances(centroids, centroids, metric='cosine')
    normal = cdist.mean()
    d2min = np.partition(dists, kth=2)[:, :2]

    a = d2min[:, 1]
    b = d2min[:, 0]
    s = a / b

    if score_method == 'median':
        score = np.median(s)
    else:
        score = np.mean(s)

    return score


def cluster_types(labels):
    """

    :param labels:
    :return:
    """
    _, counts = np.unique(labels, return_counts=True)
    thresh = 0.05 * np.mean(counts)
    counts = np.rint(counts / thresh) * thresh
    var = np.std(counts)
    ntypes = np.unique(counts).shape[0]
    return var, ntypes


def plotScores(pdb, n_range, scores, ntypes):
    """

    :param pdb:
    :param n_range:
    :param scores:
    :param ntypes:
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 10}

    # _, _, title = getPDB(pdb)
    title = pdb

    matplotlib.rc('font', **font)

    print('Plotting')
    fig, ax = plt.subplots(2, 1, figsize=(6, 3), sharex=True)
    fig.suptitle('k profile: ' + ' (' + pdb + ')')
    ax[0].scatter(n_range, scores, marker='D', s=15)
    ax[0].plot(n_range, scores)
    ax[1].plot(n_range, ntypes)
    ax[1].scatter(n_range, ntypes, marker='D', s=15)
    ax[0].axvline(x=n_range[np.argmax(scores)], label='Best Score = ' + str(n_range[np.argmax(scores)]), color='black')
    ax[1].axvline(x=n_range[np.argmax(scores)], label='Best Score', color='black')
    nc = str(n_range[np.argmax(scores)])
    ticks = ax[0].get_xticks()
    ticks = np.append(ticks, n_range[np.argmax(scores)])
    ax[1].set_xticks(ticks)
    ax[1].set_xlim([0, n_range[-1]])
    ax[1].set_xlabel('# Of Clusters')
    ax[0].set_ylabel('Quality' + '\n' + 'Score', rotation='horizontal', ha='center', va='center', labelpad=25)
    ax[1].set_ylabel('# Unique \n Clusters', rotation='horizontal', ha='center', va='center', labelpad=25)

    ax[0].tick_params(axis='y', labelsize=8)
    ax[1].tick_params(axis='y', labelsize=8)

    ax[0].grid()
    ax[1].grid()
    ax[0].legend()
    # fig.tight_layout()
    plt.show()
