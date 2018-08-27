import SBCcode.MonteCarlo.PICOcalGlobalLikelihood_v2_C as pcgl
import numpy as np
# import ipdb
import time

''' Define a new cumsum function
'''
# import numba as nb
from numba import jit, vectorize, guvectorize, uint8
import cumsum as cs

# @guvectorize([(uint8, uint8[:])], '()->(n,n,n)')
# def nb_cumsum(dim, x):
#     if dim == 0:
#         for i in range(1, x.shape[0]):
#             x[i, :, :] += x[i - 1, :, :]
#     elif dim == 1:
#         for i in range(1, x.shape[1]):
#             x[:, i, :] += x[:, i - 1, :]
#     elif dim == 2:
#         for i in range(1, x.shape[2]):
#             x[:, :, i] += x[:, :, i - 1]
#     else:
#         x = 0

# @jit(uint8(uint8[:,:,:], uint8), nopython=True)
# # @jit(nopython=True)
# def nb_cumsum(x, dim):
#     if dim == 0:
#         for i in range(1, x.shape[0]):
#             x[i, :, :] += x[i - 1, :, :]
#     elif dim == 1:
#         for i in range(1, x.shape[1]):
#             x[:, i, :] += x[:, i - 1, :]
#     elif dim == 2:
#         for i in range(1, x.shape[2]):
#             x[:, :, i] += x[:, :, i - 1]
#     else:
#         x = 0
#     return x

# @jit(uint8(uint8[:,:,:], uint8), nopython=True)


@jit(nopython=True, cache=True)
def nb_cumsum(x, dim):
    for i in range(1, x.shape[2]):
        x[:, :, i] += x[:, :, i - 1]
    return x


# testparams = np.zeros(pcgl.n_params)
# # testparams[:pcgl.n_Epts] = pcgl.Epts_initialguess.ravel()
# ipdb.set_trace()


# ll = pcgl.PICOcalLL(testparams)

# residuals = pcgl.PICOcalResiduals(testparams)

# rnd = np.array(np.random.rand(7, 12500, 20), dtype=np.uint8)
# rnd1 = np.array(np.random.rand(7, 20, 12500), dtype=np.uint8)
#

################################################################################
def test_cumsum():
    rnd = np.random.randint(0, 255, (7, 151196, 20), dtype=np.uint8)
    rnd1 = np.random.randint(0, 255, (7, 20, 151196), dtype=np.uint8)

    t0 = time.time()
    # rnd2 = np.array(rnd, order='F')
    for i in range(200):
        np.cumsum(rnd, axis=1, dtype=np.uint8)

    t1 = time.time()

    for i in range(200):
        np.cumsum(rnd1, axis=2, dtype=np.uint8)

    t2 = time.time()

    # a = rnd1[:, :, 0]
    nb_cumsum(rnd1, np.uint8(2))
    t3 = time.time()
    for i in range(200):
        nb_cumsum(rnd1, np.uint8(2))
    t4 = time.time()

    for i in range(200):
        cs.cumsum(rnd1, axis=2)

    t5 = time.time()

    print(t1 - t0)
    print(t2 - t1)
    print(t4 - t3)
    print(t5 - t4)

################################################################################


def test_chisq():
    nn = np.random.randint(0, 100, 200000)
    nu = np.random.randint(1, 100, 200000)

    t0 = time.time()

    for i in range(500):
        ch1 = pcgl.PoissonChisq(nn, nu)

    t1 = time.time()
    for i in range(500):
        ch2 = pcgl.PoissonChisq_0(nn, nu)

    t2 = time.time()
    print(t1 - t0)
    print(t2 - t1)

    print(np.amax(np.abs(ch1 - ch2)))
################################################################################


test_cumsum()
