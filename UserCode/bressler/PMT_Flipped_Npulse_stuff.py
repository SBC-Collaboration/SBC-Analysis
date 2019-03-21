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

runpath = "/bluearc/storage/SBC-17-data/20170719_0/"
events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]

V = [[],[],[],[],[],[]]
VwithNIM = [[],[],[],[],[],[]]
Vmax = []
VmaxwithNIM=[]
n=[]

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
        trace = trace - baseline
        pk_ind = scipy.signal.find_peaks(trace,5)
        pk_vals = [trace[k] for k in pk_ind[0]]
        Npeaks = len(pk_vals)
        n.append(Npeaks)
        othertrace = trac[i][1]
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
                
        for i in range(5):
            if Npeaks == i:
                for val in pk_vals:
                    V[i].append(val)

            if min(othertrace) < -30:
                for val in pk_vals:
                    VwithNIM[i].append(val)

#do the other plots
plt.figure()
vvals, bins, _= plt.hist(np.asarray(Vmax),110,color='r',histtype = 'step')
vnimvals, _, _ = plt.hist(np.asarray(VmaxwithNIM),bins=bins,color='b',histtype='step')
plt.title('RunType 902: 20170719_0')
plt.xlabel('V max')
plt.show

plt.figure()
for i in range(len(V)):
    plt.hist(np.asarray(V[i]),110,histtype='step')
plt.xlabel('V max')
plt.yscale('log')
plt.show

vnimvals = vnimvals[vvals>0]
vvals = vvals[vvals>0]

diff = vvals-vnimvals
perc = np.divide(vnimvals,vvals)
perc[np.isnan(perc)]=float('+inf')
perc=perc[perc<float('+inf')]

def functn(x,a,b):
    return scipy.stats.norm.cdf(x,a,b)

params, params_cov = scipy.optimize.curve_fit(functn,bins[:len(perc)],perc,p0=[50,1])

plt.figure()
plt.scatter(bins[:len(perc)],perc)
plt.plot(bins[:len(perc)],functn(bins[:len(perc)],params[0],params[1]),color='r')
plt.text(40,.75,"mu = "+str(params[0]),fontsize=15)
plt.text(40,.5,"sigma = "+str(params[1]),fontsize=15)
plt.xlabel('V max')
plt.ylabel('efficiency')

plt.show()
plt.figure()
plt.hist(n,max(n))
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