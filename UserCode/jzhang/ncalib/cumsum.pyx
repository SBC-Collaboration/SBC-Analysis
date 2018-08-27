"""
cumsum.pyx

"""

import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern void c_cumsum (unsigned char* array , int m, int n, int l, int axis)

@cython.boundscheck(False)
@cython.wraparound(False)
def cumsum(np.ndarray[unsigned char, ndim=3, mode="c"] input not None, int axis):
    cdef int m, n, l

    m, n, l = input.shape[0], input.shape[1], input.shape[2]

    c_cumsum (&input[0,0,0], m, n, l, axis)

    return None
