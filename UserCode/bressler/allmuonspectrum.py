#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 13:50:27 2020

@author: bressler
"""

import SBCcode as sbc
import matplotlib.pyplot as plt
import numpy as np
import runlistscatalogue as rlc
from os import listdir
from os.path import isfile,join
import gc
import pulse_integrator as pi
from gaincalc import get_gain

runs = rlc.bgOct2and3
scintillation = []
t = 0
for run in runs:
    m=get_gain("/bluearc/storage/SBC-17-data/",run)
    print(run)
    runpath = '/bluearc/storage/SBC-17-data/'+run
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    for x in events:
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,x)
        if e["fastDAQ"]["loaded"]:
            veto = e["fastDAQ"]["VetoCoinc"]
            fdt = e["fastDAQ"]["time"]
            t += fdt[-1]-fdt[0]
            if max(veto)> 0.1:
                dveto = np.diff(veto)
                imuon = list(dveto).index(max(dveto))
    
                tmuon = fdt[imuon]
                pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
                d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
                pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
                tracetimes = pmttracetime - pmtalign
                possibleCoincidenceTimes  = [tracetime for tracetime in
                                             tracetimes if abs(tracetime - tmuon)<10e-6]
                if len(possibleCoincidenceTimes)>0: 
                    #print(min(abs(possibleCoincidenceTimes-tmuon)))
                    possibledt = [abs(tracetime - tmuon) for tracetime in
                                                 tracetimes if abs(tracetime - tmuon)<10e-6]
                    #print(possibledt)
                    
                    possibleCoincidenceIndices = [list(tracetimes).index(tracetime) for tracetime in
                                                  tracetimes if abs(tracetime - tmuon)<10e-6]
                    coinci = possibleCoincidenceIndices[possibledt.index(min(possibledt))]
                    #print(coinci)
                    #for i in possibleCoincidenceIndices:
                    trace = np.fabs(e["PMTtraces"]["traces"][coinci][0]) 
                    #if ch0 saturated, stitch in low res channel:
                    if max(trace) == 128:
                        trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][coinci][1]))
                    dt = e["PMTtraces"]["dt"][coinci][0]
                    
                    #subtract baseline:
                    #Actually this gets done in pulse_integrator anyway
                    #baseline = np.mean(trace[0:50])
                    #trace -= baseline 
                                                
                    #integrate and convert to phe:
                    [phe,n,totInt,pktimes] = pi.SBC_pulse_integrator_bressler(trace,dt)
                    scintillation.append(phe/m)
                    print(phe/m)
        gc.collect()
plt.figure()
plt.hist(scintillation,int(np.ceil(np.sqrt(len(scintillation)))))
plt.xlabel('phe',fontsize=18)
plt.ylabel('count',fontsize=18)
plt.title("Rate: %f Hz"%(len(scintillation)/t))
plt.show()
