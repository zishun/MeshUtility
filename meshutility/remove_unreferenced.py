import numpy as np


__all__ = ['remove_unreferenced_vertices']


def remove_unreferenced_vertices(V, F):
    dim = F.shape[1]

    # # if not -1 in F[:, dim-1]:
    # referenced = np.zeros((V.shape[0],), '?')
    # referenced[F.ravel()] = True
    # idx = np.full((V.shape[0],), -1, 'i')
    # n = np.sum(referenced)
    # idx[referenced] = np.arange(n, dtype='i')
    # F = idx[F.ravel()].reshape((-1, dim))
    # V = V[referenced]
    
    referenced = np.zeros((V.shape[0]+1,), '?')  # add the last one for -1.
    referenced[F.ravel()] = True
    referenced[-1] = True  # suppose "-1" is always used.
    idx = np.full((V.shape[0],), -1, 'i')
    n = np.sum(referenced)
    idx[referenced[:-1]] = np.arange(n-1, dtype='i')  # map -1 to -1.

    F = idx[F.ravel()].reshape((-1, dim))
    V = V[referenced[:-1]]
    return V, F, idx
