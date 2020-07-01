#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 20:39:08 2020

@author: pi
"""

import numpy as np
cimport numpy as np
cimport cython
from cython.parallel import prange


DTYPE_uint8 = np.uint8
ctypedef np.uint8_t DTYPE_uint8_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef int count_above_cython(DTYPE_uint8_t [:] arr_view, DTYPE_uint8_t thresh) nogil:

    cdef int length, i, total
    total = 0
    length = arr_view.shape[0]

    for i in prange(length):
        if arr_view[i] >= thresh:
            total += 1

    return total


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def count_above(np.ndarray arr, DTYPE_uint8_t thresh):

    cdef DTYPE_uint8_t [:] arr_view = arr.ravel()
    cdef int total

    with nogil:
       total =  count_above_cython(arr_view, thresh)
    return total
