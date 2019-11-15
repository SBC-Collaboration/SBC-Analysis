#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 10:39:08 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
from random import randrange
import PMT_NIM_trig_efficiency as efficiency
import peakAndNimTiming

runpath = "/bluearc/storage/SBC-17-data/20170706_1/"
events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
baselines = []
V = [[],[],[],[],[],[]]
VwithNIM = [[],[],[],[],[],[]]
Vmax = []
VmaxwithNIM=[]
n=[]
NIMTimeDifferences = []
for event in events:
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)

    #d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
    #print(d.keys())
    tr = e["PMTtraces"]
    trac = tr["traces"]
    t0 = tr["t0"]
    dt = tr["dt"]
    pmttrig = e["fastDAQ"]["PMTtrig"]

    indices_of_high = []
    for i in range(len(trac)):
        # plot one random trace from each event
        j=randrange(len(trac))
        trace = np.fabs(trac[i][0])
        baseline = np.mean(trace[0:100])
        baselines.append(baseline)
        trace = trace - baseline
        pk_ind = scipy.signal.find_peaks(trace,5)
        pk_vals = [trace[k] for k in pk_ind[0]]
        Npeaks = len(pk_vals)
        n.append(Npeaks)
        othertrace = trac[i][1]
        tPMT = np.arange(len(trace))*dt[i][0]
        tNIM = np.arange(len(othertrace))*dt[i][1]
        timediff = peakAndNimTiming.findPeakNimDiff(trace,othertrace,tPMT,tNIM,False)
        if timediff:
            NIMTimeDifferences.append(timediff)
        """
        for j in range(len(diffs)):
            NIMTimeDifferences.append(diffs[j])
            """
        
        if i==j:
            plt.figure()
            plt.hold(True)
            x=np.arange(len(trace))*dt[i][0]
            plt.plot(x,trace)
            plt.scatter(pk_ind[0]*dt[i][0],pk_vals,s=50,c="r")
            plt.plot(x,othertrace,c='g')
            plt.xlabel('t (s)')
            plt.show
        
        if pk_vals:
            Vmax.append(max(pk_vals))
        if min(othertrace) < -30:
            if pk_vals:
                VmaxwithNIM.append(max(pk_vals))
                
        for i in range(len(V)):
            if Npeaks == i:
                for val in pk_vals:
                    V[i].append(val)

                if min(othertrace) < -30:
                    for val in pk_vals:
                        VwithNIM[i].append(val)

#do the other plots
"""
plt.figure()
plt.hist(NIMTimeDifferences,100,histtype='step')
plt.plot([np.mean(NIMTimeDifferences),np.mean(NIMTimeDifferences)],[0,1e4])
plt.yscale('log')
plt.xlabel('Time Delay (s), avg = '+str(np.mean(NIMTimeDifferences)))
plt.ylabel('Count')
plt.show
"""
"""
plt.figure()
for i in range(len(V)):
    plt.hist(np.asarray(V[i]),110,histtype='step')
plt.xlabel('V max')
plt.yscale('log')
plt.legend([str(i)+" peaks" for i in range(len(V))])
plt.grid()
plt.show

efficiency.NIM_efficiency_and_plot(Vmax,VmaxwithNIM,"only Vmax")
for j in range(len(V)):
    if V[j] and VwithNIM[j]:
        efficiency.NIM_efficiency_and_plot(V[j],VwithNIM[j],str(j)+" peaks")
"""

plt.figure()
nvals, bins, _ = plt.hist(n,max(n),normed=True,zorder=1)
ns = bins[:max(n)]

def pois(x,m):
    return scipy.stats.poisson.pmf(x,m)

params, params_cov = scipy.optimize.curve_fit(pois,ns,nvals)
upper = nvals + np.sqrt(np.diag(params_cov))
lower = nvals - np.sqrt(np.diag(params_cov))
print(params)
plt.plot(pois(ns,params[0]),c='r',zorder=2)
plt.fill_between(ns,lower,upper,color='r',alpha=0.5,zorder=10)
plt.xlabel('Number of Peaks found by scipy.signal.find_peaks()')
plt.show()


"""
for j in range(len(indices_of_high)):
    trace = np.fabs(trac[indices_of_high[j]][0])
    baseline = np.mean(trace[0:100])
    trace = trace - baseline
    pk_ind = scipy.signal.find_peaks(trace,5)
    pk_vals = [trace[k] for k in pk_ind[0]]
    plt.figure()
    plt.hold(True)
    x=np.arange(len(trace))
    plt.plot(x,trace)
    plt.scatter(pk_ind[0],pk_vals,s=50,c="r")
    plt.show
"""