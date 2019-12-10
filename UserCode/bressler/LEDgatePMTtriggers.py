#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 12:43:03 2019

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import matplotlib.pyplot as plt


dt = []
datadir = '/bluearc/storage/SBC-17-data'
run = '20170710_0'

for k in range(90):

    en = k
    mu = 4e7
    e = sbc.DataHandling.GetSBCEvent.GetEvent(datadir+'/'+run,en)
    cgate = e["fastDAQ"]["CAMgate"]
    dcam = np.diff(cgate)
    p1=e["fastDAQ"]["Piezo1"]
    fdt = e["fastDAQ"]["time"]
    runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
    
    camOnTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] < -0.5]
    camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]
    
    
    pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
    d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
    pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
    tracetimes = pmttracetime - pmtalign
    
    for t in tracetimes:
        if t > -0.15 and t < 0.1:
            lastCamOff = 0
            for i in range(len(camOffTimes)):
                if t > camOffTimes[i]:
                    lastCamOff = camOffTimes[i]
                elif t < camOffTimes[i]:
                    break
            dt.append(t-lastCamOff)

plt.figure()
d,b,_ = plt.hist(dt,800)
plt.show
print("bin width: "+str((b[1]-b[0])*1e6) + " microseconds")
