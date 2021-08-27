import scipy
import numpy as np
import sklearn.datasets

def print_matrix(matrix):
    for dp in matrix:
        for j in range(len(dp) - 1):
            print('{:.4f},'.format(dp[j]), end='')
        print('{:.4f}'.format(dp[-1]))

def print_array(array):
    for i in range(len(array) - 1):
        print('{:.4f},'.format(array[i]), end='')
    print('{:.4f}'.format(array[-1]))

def createWightedAdjacencyMatrix(observations, n):
     W = np.full((n, n), 0, dtype=np.float64)
     for i in range(n):
         W[i] = np.linalg.norm(observations - observations[i], axis=1) / 2
     W = np.exp(-W)
     np.fill_diagonal(W, val=0)
     return W

# Input data
X, _ = sklearn.datasets.make_blobs(n_samples=10, centers=3, n_features=2, shuffle=True, random_state=31)

# Normalized laplacian matrix:
print_matrix(scipy.sparse.csgraph.laplacian(createWightedAdjacencyMatrix(X, 10), True))
# jacobi
# eigenvalues, eigenvectors = spkmeans.jacobi(scipy.sparse.csgraph.laplacian(createWightedAdjacencyMatrix(X, 10), True).tolist())
#
# print_array(eigenvalues)
#
#
# print_matrix(np.array(eigenvectors).T) # Note: the transpose here....
