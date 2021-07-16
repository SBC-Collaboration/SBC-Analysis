#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 15:09:18 2020

@author: bressler
"""

import scipy
import numpy as np
from pulse_integrator import stitchTraces
import SBCcode as sbc
import matplotlib.pyplot as plt


def filterXeBCTrace(trace,dt):
    cutoff = 1.5e8
    nyq = 1/(2*dt)
    normal_cut = cutoff/nyq
    order = 3
    b,a = scipy.signal.butter(order, normal_cut, btype = 'low', analog = False)
    y = scipy.signal.filtfilt(b,a,trace)
    return y

def main():
    runpath = '/bluearc/storage/SBC-17-data/20170706_4'
    event = 0
    traceIndex = 40
    
    
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
    
    tr = e["PMTtraces"]
    trac = tr["traces"]
    dt = tr["dt"]
    
    trace = np.fabs(trac[traceIndex][0])
    if max(trace) == 128:
        trace = stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][traceIndex][1]))
    b = np.mean(trace[0:50])
    trace -= b
    # get the time step, assuming it's constant
    dt_tr = dt[traceIndex][0]
    tPMT = np.arange(len(trace))*dt_tr
    filteredTrace = filterXeBCTrace(trace,dt_tr)
    
    plt.figure()
    plt.plot(tPMT,trace)
    plt.plot(tPMT,filteredTrace)
    plt.show()

if  __name__ == "__main__":
    main()