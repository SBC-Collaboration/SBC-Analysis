#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 08:44:01 2021

@author: bressler
"""


import SBCcode as sbc
import numpy as np
import os
import runlistscatalogue as rlc
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('default')
import gc
import runlistscatalogue as rlc



runs = rlc.bgJuly1


T_byrun = dict.fromkeys(runs)

datadir = '/bluearc/storage/SBC-17-data/'

for run in runs:
    
    #print(run)
    runrawpath = datadir+run

    events = [evnt for evnt in os.listdir(runrawpath) if not os.path.isfile(os.path.join(runrawpath, evnt))]
    Nevents = len(events)
    T_byevent = np.zeros(Nevents)
    for i in range(Nevents):
        #if i > 2:
        #    break
        try:
            e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath, i, 'slowDAQ')
            T_mean = np.mean(e['slowDAQ']['T1'])
            T_byevent[i] = T_mean
        except Exception as x:
            print(x)
            break
        
    T_byrun[run] = np.mean(T_byevent)

plt.figure()
plt.hist([T_byrun[r] for r in runs], 10)
plt.xlabel('T1 [C]', fontsize=18)
plt.show()

print(np.mean([T_byrun[r] for r in runs]))