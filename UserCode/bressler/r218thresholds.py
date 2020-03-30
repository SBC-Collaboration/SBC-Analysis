#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 10:29:20 2019

@author: bressler
"""

from PICOcode.REFPROP.SeitzModel import SeitzModel
import numpy as np

T = np.arange(-180,-150,5)
print(T)

for t in T:
    print(t)
    m=SeitzModel(30,float(t),'argon')
    print(m.Q[0])
