#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:39:23 2019

@author: bressler
"""
import matplotlib.pyplot as plt
import SBCcode as sbc
import os
from runlistscatalogue import *
import numpy as np
from random import randrange

runs = cfJune28to30
for run in runs:
    runpath = '/bluearc/storage/SBC-17-data/%s'%run
    events = [evnt for evnt in os.listdir(runpath) if not os.path.isfile(os.path.join(runpath,evnt))]
    Nevents = len(events)
    for i in range(2):
        e=sbc.DataHandling.GetSBCEvent.GetEvent(runpath,i,'PMTtraces')
        traces = e['PMTtraces']['traces']
        tr = e["PMTtraces"]
        trac = tr["traces"]
        dt = tr["dt"]
        for i in range(len(trac)):
            trace = np.fabs(trac[i][0])
            b = np.mean(trace[0:50])
            othertrace = np.fabs(trac[i][1])
            
            dt_tr = dt[i][0]
            tPMT = np.arange(len(trace))*dt_tr
            
            if max(trace)==128:
                j = list(trace).index(128)
                multiplier = 128/othertrace[j]
                othertrace = [x*multiplier for x in othertrace]

                for i in range(len(trace)):
                    if trace[i] ==128:
                        trace[i] = othertrace[i]
                trace -= b
                if randrange(1000)==50:
                    plt.figure()
                    plt.plot(tPMT,trace)
                    plt.show
