#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 09:43:41 2024

@author: hep
"""

import numpy as np

def recondiff(d1, d2):
    for k1 in d1.keys():
        if k1 not in d2.keys():
            print('Key ' + k1 + ' only appears in the first dictionary')
            
    for k in d2.keys():
        if k not in d1.keys():
            print('Key ' + k + ' only appears in the second dictionary')
            continue
        if len(d1[k].shape) != len(d2[k].shape):
            print('Key ' + k + ' has different dimension in the two dictionaries')
            print(d1[k].shape)
            print(d2[k].shape)
            continue
        if np.any(np.array(d1[k].shape) != np.array(d2[k].shape)):
            print('Key ' + k + ' has different size in the two dictionaries')
            print(d1[k].shape)
            print(d2[k].shape)
            continue
        kdiff_ix = np.nonzero(d1[k]-d2[k])
        print('In key ' + k + ':')
        if len(kdiff_ix[0])==0:
            print('  All good!')
        else:
            for i_x in range(len(kdiff_ix[0])):
                ix = [ixx[i_x] for ixx in kdiff_ix]
                print('  ' + str(ix) + ':  ' +
                      str(d1[k][tuple(ix)]) + ' -> ' + str(d2[k][tuple(ix)]))
        
    
