#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 14:43:06 2020

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import matplotlib.pyplot as plt
from PICOcode.REFPROP.SeitzModel import SeitzModel
import runlistscatalogue as rlc

runs=["20170920_1"]
LT = []
setpoints = []
LTpreq=5
indices_back=30
triggerp = []
for run in runs:
    runrawpath = '/bluearc/storage/SBC-17-data/'+run
    runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
    historyfilename = runreconpath+"HistoryAnalysis_%s.bin"%run
    history = sbc.DataHandling.ReadBinary.ReadBlock(historyfilename)
    for eventn in range(10):
        try:     
            e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,eventn,"slowDAQ","event")
            T = np.mean(e["slowDAQ"]["T1"])
            t = e["slowDAQ"]["elapsed_time"]
            print(T)
            pset = e["event"]["Pset"]
            seitz = SeitzModel(pset,T,'xenon')
            Q = seitz.Q
            PT6 = e["slowDAQ"]["PT6"]
            PT5 = e["slowDAQ"]["PT5"]
            plt.figure()
            PT4 = e["slowDAQ"]["PT4"]
            plt.plot(t,PT6,label="PT6")
            plt.plot(t,PT5,label="PT5")
            plt.plot(t,PT4,label="PT4")
            plt.plot([t[0],t[-1]],[pset,pset],label="pset")
            plt.legend()
            print("event livetime: "+str(e["event"]["livetime"]))
            trig_pressure = e["slowDAQ"]["PT6"][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
            
            edges = history["PressureEdge"][5]
            times = history["PressureBins"][eventn][5][:]
            centers = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
            if pset not in setpoints:
                setpoints.append(pset)
                LT.append(0)
                triggerp.append([])
            ind = setpoints.index(pset)
            triggerp[ind].append(abs(trig_pressure-pset))
            s = 0
            for i in range(len(centers)):
                if abs(centers[i]-pset)<LTpreq:
                    s += times[i]
                    LT[ind]+=times[i]
            print("counted time: "+str(s))

        except Exception as x:
            print(x)
            break
for i in range(len(triggerp)):
    plt.figure()
    plt.hist(triggerp[i],20)
    plt.title(str(setpoints[i]))
    plt.show