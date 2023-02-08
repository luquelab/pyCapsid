"""Module with functions for identifying quasi-rigid clusters given a distance fluctuation matrix, either calculated using
ENM or provided."""
import numpy as np

def findQuasiRigidClusters(pdb, dist_flucts, n_range, cluster_method='discretize', return_type='final', score_method='median', save=False, dir='.'):
    """

    :param str pdb:
    :param dist_flucts:
    :param list n_range:
    :param str cluster_method:
    :param str score_method:
    :param bool save:
    :param str dir:
    :return:
    """
    sims = fluctToSims(dist_flucts)
    n_vecs = n_range.max()

    embedding = calcEmbedding(sims, n_vecs)

    labels, scores, variances, numtypes = cluster_embedding(n_range, embedding, method=cluster_method, score_method=score_method)

    from .clustering_util import plotScores
    plotScores(pdb, n_range, scores, numtypes)

    ind = np.argmax(scores)
    final_clusters = labels[ind]
    final_score = scores[ind]
    if return_type=='final':
        return final_clusters, final_score
    elif return_type=='full':
        return labels, scores, variances, numtypes
    else:
        return final_clusters

def fluctToSims(d):
    """Transforms a distance fluctuation matrix into a similarity matrix

    :param d: Sparse distance fluctuation matrix
    :return: Sparse similarity matrix
    """
    from scipy import sparse
    d_bar = np.mean(np.sqrt(d.data))
    # print('RMS distance fluctuations: ', d_bar)
    sigma = 1 / (2 * d_bar ** 2)
    data = d.data
    data = np.exp(-sigma * data ** 2)
    sims = sparse.coo_matrix((data, (d.row, d.col)), shape=d.shape)
    sims.eliminate_zeros()
    return sims


def calcEmbedding(sims, n_vecs):
    """

    :param sims:
    :param n_vecs:
    :return:
    """
    from sklearn.manifold import spectral_embedding
    print('Performing Spectral Embedding')


    X_transformed = spectral_embedding(sims, n_components=n_vecs, drop_first=False, eigen_solver='arpack',
                                       norm_laplacian=True)

    return X_transformed


def cluster_embedding(n_range, maps, method='discretize', score_method='median'):
    """

    :param n_range:
    :param maps:
    :param method:
    :param score_method:
    :return:
    """
    print('Clustering Embedded Points')
    print(f'Method: {method}')

    from sklearn.preprocessing import normalize
    if method == 'kmeans':
        from sklearn.cluster import k_means
    # from sklearn.cluster import MiniBatchKMeans, AgglomerativeClustering
    # from sklearn.metrics import pairwise_distances
    # from sklearn_extra.cluster import KMedoids

    # from sklearn.metrics import silhouette_score
    # from sklearn.metrics import davies_bouldin_score
    from .clustering_util import median_score, cluster_types, calcCentroids, calcCosCentroids

    labels = []
    scores = []
    variances = []
    numtypes = []

    for n in range(len(n_range)):
        n_clusters = n_range[n]
        emb = maps[:, :n_clusters].copy()

        normalize(emb, copy=False)

        # print('Clusters: ' + str(n_clusters))

        if method == 'discretize':
            label = discretize(emb)
            centroids = calcCosCentroids(emb, label, n_clusters)
        elif method == 'kmeans':
            centroids, label, inert, n_iter = k_means(emb, n_clusters=n_clusters,
                                                      return_n_iter=True)
        else:
            raise Exception('Method should be kmeans or discretize.')

        labels.append(label)

        testScore = median_score(emb, centroids, score_method)
        # print('Score: ', testScore)
        scores.append(testScore)

        var, ntypes = cluster_types(label)
        variances.append(var)
        numtypes.append(ntypes)

    return labels, scores, variances, numtypes


# Adapted from sklearn. May need to change later.
def discretize(
        vectors, *, copy=False, max_svd_restarts=300, n_iter_max=2000, random_state=None
):
    """Search for a partition matrix which is closest to the eigenvector embedding.
    This implementation was proposed in [1]_.
 
    :param vectors: array-like of shape (n_samples, n_clusters)
        The embedding space of the samples.
    :param copy: bool, default=True
        Whether to copy vectors, or perform in-place normalization.
    :param max_svd_restarts: int, default=30
        Maximum number of attempts to restart SVD if convergence fails
    :param n_iter_max: int, default=30
        Maximum number of iterations to attempt in rotation and partition
        matrix search if machine precision convergence is not reached
    :param random_state: int, RandomState instance, default=None
        Determines random number generation for rotation matrix initialization.
        Use an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.
    Returns
    -------
    labels: array of integers, shape: n_samples
        The labels of the clusters.
    References
    ----------
    .. [1] `Multiclass spectral clustering, 2003
           Stella X. Yu, Jianbo Shi
           <https://www1.icsi.berkeley.edu/~stellayu/publication/doc/2003kwayICCV.pdf>`_
    Notes
    -----
    The eigenvector embedding is used to iteratively search for the
    closest discrete partition.  First, the eigenvector embedding is
    normalized to the space of partition matrices. An optimal discrete
    partition matrix closest to this normalized embedding multiplied by
    an initial rotation is calculated.  Fixing this discrete partition
    matrix, an optimal rotation matrix is calculated.  These two
    calculations are performed until convergence.  The discrete partition
    matrix is returned as the clustering solution.  Used in spectral
    clustering, this method tends to be faster and more robust to random
    initialization than k-means.
    """

    from scipy.sparse import csc_matrix
    from scipy.linalg import LinAlgError
    from sklearn.utils import check_random_state, as_float_array

    random_state = check_random_state(random_state)

    vectors = as_float_array(vectors, copy=copy)

    eps = np.finfo(float).eps
    n_samples, n_components = vectors.shape

    # Normalize the eigenvectors to an equal length of a vector of ones.
    # Reorient the eigenvectors to point in the negative direction with respect
    # to the first element.  This may have to do with constraining the
    # eigenvectors to lie in a specific quadrant to make the discretization
    # search easier.
    norm_ones = np.sqrt(n_samples)
    for i in range(vectors.shape[1]):
        vectors[:, i] = (vectors[:, i] / np.linalg.norm(vectors[:, i])) * norm_ones
        if vectors[0, i] != 0:
            vectors[:, i] = -1 * vectors[:, i] * np.sign(vectors[0, i])

    # Normalize the rows of the eigenvectors.  Samples should lie on the unit
    # hypersphere centered at the origin.  This transforms the samples in the
    # embedding space to the space of partition matrices.
    vectors = vectors / np.sqrt((vectors ** 2).sum(axis=1))[:, np.newaxis]

    svd_restarts = 0
    has_converged = False

    # If there is an exception we try to randomize and rerun SVD again
    # do this max_svd_restarts times.
    while (svd_restarts < max_svd_restarts) and not has_converged:

        # Initialize first column of rotation matrix with a row of the
        # eigenvectors
        rotation = np.zeros((n_components, n_components))
        rotation[:, 0] = vectors[random_state.randint(n_samples), :].T

        # To initialize the rest of the rotation matrix, find the rows
        # of the eigenvectors that are as orthogonal to each other as
        # possible
        c = np.zeros(n_samples)
        for j in range(1, n_components):
            # Accumulate c to ensure row is as orthogonal as possible to
            # previous picks as well as current one
            c += np.abs(np.dot(vectors, rotation[:, j - 1]))
            rotation[:, j] = vectors[c.argmin(), :].T

        last_objective_value = 0.0
        n_iter = 0

        while not has_converged:
            n_iter += 1

            t_discrete = np.dot(vectors, rotation)

            labels = t_discrete.argmax(axis=1)
            vectors_discrete = csc_matrix(
                (np.ones(len(labels)), (np.arange(0, n_samples), labels)),
                shape=(n_samples, n_components),
            )

            t_svd = vectors_discrete.T * vectors

            try:
                U, S, Vh = np.linalg.svd(t_svd)
            except LinAlgError:
                svd_restarts += 1
                print("SVD did not converge, randomizing and trying again")
                break

            ncut_value = 2.0 * (n_samples - S.sum())
            if (abs(ncut_value - last_objective_value) < eps) or (n_iter > n_iter_max):
                has_converged = True
            else:
                # otherwise calculate rotation and continue
                last_objective_value = ncut_value
                rotation = np.dot(Vh.T, U.T)

    if not has_converged:
        raise LinAlgError("SVD did not converge")
    return labels


# Also adapted from sklearn, similar method???
# def cluster_qr(vectors):
#     from scipy.linalg import svd, qr
#     """Find the discrete partition closest to the eigenvector embedding.
#         This implementation was proposed in [1]_.
#     .. versionadded:: 1.1
#         Parameters
#         ----------
#         vectors : array-like, shape: (n_samples, n_clusters)
#             The embedding space of the samples.
#         Returns
#         -------
#         labels : array of integers, shape: n_samples
#             The cluster labels of vectors.
#         References
#         ----------
#         .. [1] :doi:`Simple, direct, and efficient multi-way spectral clustering, 2019
#             Anil Damle, Victor Minden, Lexing Ying
#             <10.1093/imaiai/iay008>`
#     """
#
#     k = vectors.shape[1]
#     _, _, piv = qr(vectors.T, pivoting=True)
#     ut, _, v = svd(vectors[piv[:k], :].T)
#     vectors = abs(np.dot(vectors, np.dot(ut, v.conj())))
#     return vectors.argmax(axis=1)
