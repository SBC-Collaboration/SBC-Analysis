#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 08:01:25 2021

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


run1 = '20170919_0'
run2 = '20170924_0'
run3 = '20170926_3'
run4 = '20171001_0'
run5 = '20171002_0'
run6 = '20171006_0'
run7 = '20171009_8'
run8 = '20170923_3'
run9 = '20170927_4'
run10 = '20170930_1'
run11 = '20170926_4'
run12 = '20170928_0'
run13 = '20170928_5'
run14 = '20170930_7'
run15 = '20171003_4'
run16 = '20171004_0'
run17 = '20171007_0'
run18 = '20171009_0'
run19 = '20171010_7'
run20 = '20171011_8'

augrun1 = '20170801_0'
augrun2 = '20170802_0'
augrun3 = '20170803_0'
augrun4 = '20170804_0'
augrun5 = '20170805_0'


runs = [run1, run2, run3, run4, run5, run6, run7, run8, run9, run10, run11, 
        run12, run13, run14, run15, run16, run17, run18, run19, run20]

augruns = [augrun1, augrun2, augrun3, augrun4, augrun5]


Pdiff = dict.fromkeys(runs)
T5 = dict.fromkeys(runs)

Pdiff_fill_sept = -0.8
T5_fill_sept = -98.57

Pdiff_fill_aug = 7.5
T5_fill_aug = -95.57

Pdiff_aug = dict.fromkeys(augruns)
T5_aug = dict.fromkeys(augruns)
datadir = '/bluearc/storage/SBC-17-data/'

for run in runs:
    
    print(run)
    runrawpath = datadir+run

    events = [evnt for evnt in os.listdir(runrawpath) if not os.path.isfile(os.path.join(runrawpath, evnt))]
    Nevents = len(events)
    differences = []
    T5_byevent = np.zeros(Nevents)
    for i in range(Nevents):
        #if i > 2:
        #    break
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath, i, 'slowDAQ')
        
        PT4 = e['slowDAQ']['PT4']
        PT6 = e['slowDAQ']['PT6']
        for k in range(len(PT4)):
            if abs(PT4[k]-e['slowDAQ']['PressureSetpoint'][k]) < 2:
                differences.append((PT4[k]-PT6[k])-Pdiff_fill_sept)
                
        T5_mean = np.mean(e['slowDAQ']['T5']) - T5_fill_sept
        T5_byevent[i] = T5_mean
    Pdiff[run] = differences
    T5[run] = np.mean(T5_byevent)
    
    
for run in augruns:
    
    print(run)
    runrawpath = datadir+run

    events = [evnt for evnt in os.listdir(runrawpath) if not os.path.isfile(os.path.join(runrawpath, evnt))]
    Nevents = len(events)
    differences = []
    for i in range(Nevents):
        #if i > 2:
        #    break
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath, i, 'slowDAQ')
        
        PT4 = e['slowDAQ']['PT4']
        PT6 = e['slowDAQ']['PT6']
        for k in range(len(PT4)):
            if abs(PT4[k]-e['slowDAQ']['PressureSetpoint'][k]) < 2:
                differences.append((PT4[k]-PT6[k])-Pdiff_fill_aug)
                
        T5_mean = np.mean(e['slowDAQ']['T5']) - T5_fill_aug

    Pdiff_aug[run] = differences
    T5_aug[run] = T5_mean
    
plt.figure()
for run in runs:
    print(np.mean(Pdiff[run]))
    print(T5[run])
    print(len(Pdiff[run]))
    plt.hist(Pdiff[run], int(np.ceil(abs(np.sqrt(len(Pdiff[run]))))), density=True, linewidth=3,  histtype='step', label=run)
plt.xlabel('PT4 - PT6', fontsize=18)
plt.ylabel('density')
plt.legend(fontsize=18)
plt.show()

T5_points = [T5[run] for run in runs]
Pdiff_points = [np.mean(Pdiff[run]) for run in runs]
T5_points_aug = [T5_aug[run] for run in augruns]
Pdiff_points_aug =[np.mean(Pdiff_aug[run]) for run in augruns]

p, pcov = curve_fit(lambda x,m,b : m*x + b, T5_points, Pdiff_points, [-1, 0])

plt.figure()
plt.scatter(T5_points, Pdiff_points, 30, 'b', label='September fill')
plt.scatter(T5_points_aug, Pdiff_points_aug, 30, 'k', label='August fill')
plt.plot(np.arange(-10,0,0.1), [p[0]*x + p[1] for x in np.arange(-10,0,0.1)], 'b', label=r'y=%.2f x + %.2f'%(p[0], p[1]))
plt.xlabel(r'$T_{run}-T_{fill}$ [Celsius]', fontsize=20)
plt.ylabel(r'$\Delta P_{run} - (PT8-PT6)_{fill}$', fontsize=20)
#plt.ylabel(r'$\Delta P_{run} - \Delta P_{fill}$', fontsize=20)

plt.legend(fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()