#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 14:27:34 2019

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import matplotlib.pyplot as plt
from pulse_integrator import SBC_pulse_integrator_bressler


def PMTandPiezoPlot(datadir,run,event,gain):
    """
    Plots the piezo trace from piezo 0 and PMT pulse sizes near the acoustic
        signal on one figure. The green line indicates the time of acoustic t0.
        The red line indicates the time and size of the PMT t0 chosen by Matt's
        analysis. The yellow lines indicate the times and sizes of other PMT
        pulses.
    Also plots in separate windows the PMT traces for pulses 500 microseconds
        or less before acoustic t0.
    """
    en = event
    mu = gain
    e = sbc.DataHandling.GetSBCEvent.GetEvent(datadir+'/'+run,en)
    p1=e["fastDAQ"]["Piezo1"]
    fdt = e["fastDAQ"]["time"]
    runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
    pmtdiffs = []
    diffs = []
    
    acousticfilename = runreconpath+"AcousticAnalysis_%s.bin"%run
    a = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
    bubt0 = a["bubble_t0"]
    
    pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
    d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
    pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
    tracetimes = pmttracetime - pmtalign
    at0 = bubt0[en,0]
    i=0
    candidates = []
    candidate_times=[]
    pmtt0=0
    for t in (tracetimes-at0):
    
        if t<0.2 and t>-0.2:
    
            pmtdiffs.append(t)
            trace = -(e["PMTtraces"]["traces"][i][0])
            dt = e["PMTtraces"]["dt"][i][0]
            baseline = np.mean(trace[0:50])
            trace = trace - baseline
            tPMT = np.arange(len(trace))*dt
            [phe,n,totInt,pktimes] = SBC_pulse_integrator_bressler(trace,dt)
            
            if phe != None:
                phe /= mu
                candidates.append(phe)
                candidate_times.append(t)
            i+=1
    candidate_phe = 0
    the_index = 0
    i=0
    near_trace_indices = []
    for t in candidate_times:
        if t > -500e-6 and t <0:
            near_trace_indices.append(i)
            if candidates[i]>candidate_phe:
                candidate_phe = candidates[i]
                the_index = i
        i+=1
            
    if len(candidates) != 0:
        if max(candidates)>0:
            diffs.append(candidate_times[candidates.index(max(candidates))])
    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(fdt,p1,'b')
    for i in range(len(candidates)):
        if i == the_index:
            ax2.plot([candidate_times[i]+at0,candidate_times[i]+at0],[0,candidates[i]],'r',lw=4)
        else:
            ax2.plot([candidate_times[i]+at0,candidate_times[i]+at0],[0,candidates[i]],'y',lw=4)
    ax2.plot([min(candidate_times),max(candidate_times)],[0,0],linewidth=2)
    ax2.plot([at0,at0],[0,max(candidates)],'g',linewidth=4)
    ax2.set_ylabel('pmt signal (phe)',fontsize=20)
    ax1.set_xlabel('time (s)',fontsize=20)
    ax1.set_ylabel('Acoustic signa(V)',fontsize=20)
    ax1.set_ylim([min(p1),max(p1)])
    ax2.set_xlim([-0.03,-0.015])
    plt.show
    
    for i in near_trace_indices:
        trace = e["PMTtraces"]["traces"][i][0]
        dt = e["PMTtraces"]["dt"]
        dt_tr = dt[i][0]
        tPMT = np.arange(len(trace))*dt_tr
        plt.figure()
        plt.plot(tPMT,trace)
        plt.xlabel("t (s)")
        plt.ylabel("PMT signal")
        plt.show
    """
    plt.figure()
    plt.plot(e["fastDAQ"]["time"],e["fastDAQ"]["VetoCoinc"])
    plt.show
    """
def main():
    PMTandPiezoPlot('/bluearc/storage/SBC-17-data', "20170710_4",42, 4e7)
    
if __name__=="__main__":
    main()