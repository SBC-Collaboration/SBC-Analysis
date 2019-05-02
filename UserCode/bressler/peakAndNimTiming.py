#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:04:47 2019

@author: bressler
"""
import scipy
import numpy as np
import matplotlib.pyplot as plt
import SBCcode as sbc


def findPeakNimDiff(PMTtrace,NIMtrace,timePMT,timeNIM,plotbool):
    baseline = np.mean(PMTtrace[0:100])
    PMTtrace = PMTtrace - baseline
    pk_ind = scipy.signal.find_peaks(PMTtrace,5)
    pk_vals = [PMTtrace[k] for k in pk_ind[0]]
    pk_times = [timePMT[k] for k in pk_ind[0]]
    Npeaks = len(pk_vals)
    
    dt_NIM = timeNIM[1]-timeNIM[0]
    pos_derivative = -np.diff(NIMtrace)/dt_NIM
    NIMDerivativePeaks = scipy.signal.find_peaks(pos_derivative,1e10)
    #print(len(NIMDerivativePeaks[0]))
    
    if len(NIMDerivativePeaks[0])>0:
        NIMPulseTime = NIMDerivativePeaks[0][0]*dt_NIM
    else: 
        NIMPulseTime = False 
    if plotbool:
        plt.figure()
        plt.plot(timePMT,PMTtrace)
        plt.plot(timeNIM,NIMtrace)
        #plt.plot(timeNIM[0:len(pos_derivative)],pos_derivative)
        plt.scatter(pk_times,pk_vals,c='r',s=50)
        plt.scatter(NIMPulseTime,0,c='r',s=100)
        plt.show
    
    if Npeaks == 1 and NIMPulseTime:
        if pk_times:
            return pk_times[0]-NIMPulseTime
    
    #return [NIMPulseTime-pk_times[i] for i in range(len(pk_times))]
    
def main():
    e = sbc.DataHandling.GetSBCEvent.GetEvent("/bluearc/storage/SBC-17-data/20170719_0/",0)

    #d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
    #print(d.keys())
    tr = e["PMTtraces"]
    trac = tr["traces"]
    dt = tr["dt"]

    trace = np.fabs(trac[0][0])
    baseline = np.mean(trace[0:100])
    trace = trace - baseline
    othertrace = trac[0][1]
    tPMT = np.arange(len(trace))*dt[0][0]
    tNIM = np.arange(len(othertrace))*dt[0][1]
    timediff = findPeakNimDiff(trace,othertrace,tPMT,tNIM,True)
    
if __name__ == "__main__":
    main()
    