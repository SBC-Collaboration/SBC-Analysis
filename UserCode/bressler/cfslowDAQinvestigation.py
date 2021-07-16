#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 08:54:40 2020

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import scipy
import matplotlib.pyplot as plt
from collections import Counter
import runlistscatalogue as rlc
import matplotlib.pyplot as plt
from PICOcode.REFPROP.SeitzModel import SeitzModel
from plot_xebc_slowDAQ import plot_xebc_slowDAQ

runs = rlc.CfCombinedMultiTemp
Qbins=[1,1.5,2,2.5,3]
evids = [[],[],[],[]]
counts = []
nbubs = []
bubInfo = []
eventcount = np.zeros(87)
spectra = [[] for i in range(87)]
LT = np.zeros(87)
LT_by_Q = np.zeros(len(Qbins)-1)
expand_times = [[] for i in range(87)]
totbub = 0
onebubcount = np.zeros(87)
twobubcount = np.zeros(87)
threebubcount = np.zeros(87)
setpoints = []
elt = 0
temps = []
thresholds = []
pressures = []
triggers = []
for run in runs:
    print(run)
    tcut = 20
    z_low = -3
    z_high = 0
    preq = 1
    runrawpath = '/bluearc/storage/SBC-17-data/'+run
    runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
    historyfilename = runreconpath+"HistoryAnalysis_%s.bin"%run
    history = sbc.DataHandling.ReadBinary.ReadBlock(historyfilename)
    getbubfile = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/HumanGetBub_%s.bin"%run
    c = sbc.DataHandling.ReadBinary.ReadBlock(getbubfile)
    eventn = c["ev"]
    count = Counter(eventn)
    edges = history["PressureEdge"][5]
    centersp = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
    #centers = Qvals
    e0 = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,0,"slowDAQ")
    T=np.mean(e0["slowDAQ"]["T1"])
    
    sm = SeitzModel(list(edges),T,'xenon')
    edges_by_Q = sm.Q
    centers = [(edges_by_Q[i] + edges_by_Q[i+1])/2 for i in range(len(edges_by_Q)-1)]

    print("T=%f C"%T)
    
    for eventn in range(101):
        try:
            indices_back = 30
            event_ID = run+"-"+str(eventn)
            e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,eventn,"slowDAQ","event")
            sd = e["slowDAQ"]
            t =sd["elapsed_time"]
            T = np.mean(sd["T1"])
            pset = e["event"]["Pset"]
            tcenters = [(t[i+1] + t[i])/2 for i in range(len(t)-1)]
            trig_time = t[list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
            pslope = np.diff(e["slowDAQ"]["PT6"])
            expstartind = list(pslope).index(min(list(pslope)))
            exp_start = tcenters[expstartind]
            trig_pressure = e["slowDAQ"]["PT6"][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
            
            SM = SeitzModel(trig_pressure,T,'xenon')
            Q = (SM.Q)[0]
            if (not np.isnan(Q)) and (trig_time > exp_start + tcut) and(abs(trig_pressure-pset)<preq):
                temps.append(T)
                thresholds.append(Q)
                pressures.append(trig_pressure)
                for j in range(len(Qbins)-1):
                    if Q>Qbins[j] and Q<Qbins[j+1]: evids[j].append(event_ID)
                
        except Exception as x:
            print(x)
print(len(evids[1]))
plt.figure()
plt.hist(thresholds,int(np.ceil(np.sqrt(len(thresholds)))))
plt.show()

plt.figure()
plt.hist(temps,int(np.ceil(np.sqrt(len(temps)))))
plt.show()

plt.figure()
plt.hist(pressures,int(np.ceil(np.sqrt(len(pressures)))))
plt.show()
"""
for evid in evids[1]:
    runpath = '/bluearc/storage/SBC-17-data/'+evid.split('-')[0]
    n = int(evid.split('-')[1])
    plot_xebc_slowDAQ(runpath,n)
"""
