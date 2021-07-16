#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 09:20:45 2020

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import runlistscatalogue as rlc
from PICOcode.REFPROP.SeitzModel import SeitzModel
import pulse_integrator as pi
import gc

runs = rlc.cfJune28to30
Qbins = [1,1.5,2,2.5,3,3.5]


eventcount = np.zeros(87)
spectra = [[] for i in range(87)]
LT = np.zeros(87)
LT_by_Q = np.zeros(len(Qbins)-1)
expand_times = [[] for i in range(87)]
onebubcount = np.zeros(87)
twobubcount = np.zeros(87)
threebubcount = np.zeros(87)
Qvals = [[] for i in range(87)]
setpoints = []
elt = 0
temps = []
m=4e7
for run in runs:
    print(run)
    tcut = 0
    
    runrawpath = '/bluearc/storage/SBC-17-data/'+run
    runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
    historyfilename = runreconpath+"HistoryAnalysis_%s.bin"%run
    history = sbc.DataHandling.ReadBinary.ReadBlock(historyfilename)
    
    edges = history["PressureEdge"][5]
    centersp = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
    #centers = Qvals
    e0 = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,0,"slowDAQ")
    T=np.mean(e0["slowDAQ"]["T1"])
    temps.append(T)
    sm = SeitzModel(list(edges),T,'xenon')
    edges_by_Q = sm.Q
    centers = [(edges_by_Q[i] + edges_by_Q[i+1])/2 for i in range(len(edges_by_Q)-1)]

    print("T=%f C"%T)

    
    for eventn in range(101):
        try:
            indices_back = 30
            event_ID = run+"-"+str(eventn)
            e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,eventn)
            acousticfilename = runreconpath+"AcousticAnalysis_%s.bin"%run
            a = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
            bubt0 = a["bubble_t0"]

            t = e["slowDAQ"]["elapsed_time"]
            tcenters = [(t[i+1] + t[i])/2 for i in range(len(t)-1)]
            trig_time = t[list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
            pslope = np.diff(e["slowDAQ"]["PT6"])
            expstartind = list(pslope).index(min(list(pslope)))
            exp_start = tcenters[expstartind]
            trig_pressure = e["slowDAQ"]["PT6"][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
            
            SM = SeitzModel(trig_pressure,T,'xenon')
            Q = (SM.Q)[0]
            #print(Q[0])
            pset = e["event"]["Pset"]
            ev_lt = e["event"]['livetime']
            elt += ev_lt
            if pset not in setpoints:
                setpoints.append(pset)
            #print(len(history["PressureBins"][eventn]))
            times = history["PressureBins"][eventn][5][:]
            for i in range(len(centers)):
                if trig_time > exp_start + tcut:
                    if trig_pressure >= edges[i] and trig_pressure <= edges[i+1]:
                        theIndexIWant = i
                    #LT[i] += times[i]
                    #for j in range(len(Qbins)-1):
                        #if Q >= Qbins[j] and Q <= Qbins[j+1]:
                            #LT_by_Q[j] += times[i]
            
            cgate = e["fastDAQ"]["CAMgate"]
            dcam = np.diff(cgate)
            fdt = e["fastDAQ"]["time"]
            camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]

            pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
            d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
            pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
            tracetimes = pmttracetime - pmtalign
            at0 = bubt0[int(eventn),0]
            i=1 # to match the indexing of the pre-made code 
            candidate = 0
            candidate_time = 0
            candidate_PMTtime = 0
            candidate_index = 0
            active=0.01
            for t in (tracetimes-at0):
                if t<0 and t>-active:
                    trace = np.fabs(e["PMTtraces"]["traces"][i][0]) 
                    #if ch0 saturated, stitch in low res channel:
                    if max(trace) == 128:
                        trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][i][1]))
                        dt = e["PMTtraces"]["dt"][i][0]
                    for i in range(len(centers)):
                        if trig_time > exp_start + tcut:
                            if trig_pressure >= edges[i] and trig_pressure <= edges[i+1]:
                                LT[i] += active
                            for j in range(len(Qbins)-1):
                                if Q >= Qbins[j] and Q <= Qbins[j+1]:
                                    LT_by_Q[j] += active
                                            
                    #integrate and convert to phe:
                    [a,n,totInt,pktimes] = pi.SBC_pulse_integrator_bressler(trace,dt)
                    light = a/m
                    spectra[theIndexIWant].append(light)
                    Qvals[theIndexIWant].append(Q)
                    gc.collect()
                i += 1
        except Exception as x:
            print(x)
hy = []
hx = []
for i in range(len(expand_times)):
    if len(spectra[i])>0:
        for j in range(len(spectra[i])):
            hx.append(centersp[i])
            hy.append(spectra[i][j])
plt.figure()
plt.hist2d(hx,hy,bins=(40,int(max(hy))),cmap='Greys')
plt.colorbar()
plt.xlabel('PT6 [psia]',fontsize=18)
plt.ylabel('Light collected [phe]',fontsize=18)
plt.yscale('symlog',linthreshy = 0.9)
plt.show


xedges = np.arange(88)
#yedges = np.arange(int(max(hy))+1)
bins = [2**i+0.5 for i in range(11)]
#bins = np.arange(0.5,1+np.ceil(max(spect)))
bins = np.insert(bins,0,0.5)
bins=np.insert(bins,0,0)
bins=np.insert(bins,0,-100)
H, xedges, yedges = np.histogram2d(hx,hy,bins=(xedges,bins))
for i in range(len(xedges)-1):
    #print(hy[i])
    for j in range(len(H[i,:])):
        if LT[i] > 0:
            H[i,j] = H[i,j]/LT[i]
        else:
            H[i,j] = 0
print("max rate: %f"%np.amax(H))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
X,Y = np.meshgrid(xedges,yedges)

im=ax.pcolormesh(X,Y,H.T,cmap="magma")
plt.yscale('symlog',linthreshy = 0.9)
fig.colorbar(im,ax=ax)
plt.show()

hy = []
hx = []
for i in range(len(expand_times)):
    if len(spectra[i])>0:
        for j in range(len(spectra[i])):
            hx.append(Qvals[i][j])
            hy.append(spectra[i][j])

Qedges = Qbins    
#print(Qbins)
#yedges = np.arange(int(max(hy))+1)
bins = [(2**i)+0.5 for i in range(11)]
#bins = np.arange(0.5,1+np.ceil(max(spect)))
bins = np.insert(bins,0,0.5)
bins=np.insert(bins,0,-0.5)
bins=np.insert(bins,0,-100)
binc=[(bins[i+1]+bins[i])/2 for i in range(len(bins)-1)]
H2, Qedges, yedges = np.histogram2d(hx,hy,bins=(Qedges,bins))
ns = H2.copy()
errorbararray = np.zeros_like(H2)
print("total number of events in H2: %d"%np.sum(H2))
for i in range(len(LT_by_Q)):
    for j in range(len(H[i,:])):
        print(LT_by_Q[i])
        H2[i,j] = H2[i,j]/LT_by_Q[i]
        errorbararray[i,j]= np.sqrt(ns[i,j])/LT_by_Q[i]
        #print(H[i,j])
        #print(errorbararray[i,j])
        if np.isnan(H2[i,j]):
            H2[i,j]=0
            errorbararray[i,j] = 0
print("Max in H2: %f"%np.amax(H2))
print("total scintillation rate: %f"%np.sum(H2))
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.title.set_text('Rates [Hz/bin]')
X,Y = np.meshgrid(Qedges,yedges)

im=ax.pcolormesh(X,Y,H2.T,cmap="Greys")
plt.yscale('symlog',linthreshy = 0.9)
plt.ylabel("Collected light [phe]")
plt.xlabel("Seitz Threshold [keV]")
fig.colorbar(im,ax=ax)

ax2 = fig.add_subplot(2,1,2)
ax2.title.set_text('Number of bubbles in each bin')
im2=ax2.pcolormesh(X,Y,ns.T,cmap="Blues")
plt.yscale('symlog',linthreshy = 0.9)
plt.ylabel("Collected light [phe]")
plt.xlabel("Seitz Threshold [keV]")
fig.colorbar(im2,ax=ax2)
#plt.show()

fig=plt.figure()
ax3 = fig.add_subplot(1,1,1)
ax3.title.set_text("Rates per Scintillation")
for k in range(len(Qedges)-1):
    ax3.errorbar(binc,H2[k,:],errorbararray[k,:],marker='.',markersize=7,
                 linestyle='none',label="%f-%f keV"%(Qedges[k],Qedges[k+1]))
plt.xscale('symlog',linthreshx=0.5)
plt.yscale('symlog',linthreshy=0.001)
plt.legend(fontsize=12)
plt.xlabel("collected light [phe]",fontsize=18)
plt.ylabel("Rate [Hz]",fontsize=18)
#plt.ylim([0,1])
plt.grid(which='both',axis='both')
plt.show()
