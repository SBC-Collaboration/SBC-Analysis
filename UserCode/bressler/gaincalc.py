#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:00:06 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
from pulse_integrator import SBC_pulse_integrator_bressler

def get_gain(datapath,run):
    
    runpath = datapath+"/"+run+'/'
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    allTraces = []
    notrig_pulses = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    for event in events:
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)

        tr = e["PMTtraces"]
        trac = tr["traces"]
        dt = tr["dt"]
        for i in range(len(trac)):
            trace = np.fabs(trac[i][0])
            b = np.mean(trace[0:100])
    
            # get the time step, assuming it's constant
            dt_tr = dt[i][0]

            # populate dictionaries arrays based on how many pulses there were

            [a,n,totInt,pktimes] = SBC_pulse_integrator_bressler(trace,dt_tr)
            if n == 0:
                number = 'zero'
            elif n == 1:
                number = 'one'
            elif n == 2:
                number = 'two'
            elif n == 3:
                number = 'three'
            else:
                number = 'other'
            notrig_pulses[number].append(a)
            allTraces.append(a)

    for k in notrig_pulses:
        notrig_pulses[k] = [x for x in notrig_pulses[k] if x != None]
    
    allTraces = [x for x in allTraces if x != None]
    
    plt.figure()
    Nbins = int(np.floor(np.sqrt(len(allTraces))))
    allvals, bins, _ = plt.hist(allTraces,Nbins,label='all traces')
    areaVals_notrig = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    
    for k in notrig_pulses:
        areaVals_notrig[k], _, _ = plt.hist(notrig_pulses[k],bins,histtype = 'step',linewidth = 3,label='no trig '+k+' hits')
        
    spe_spectrum = areaVals_notrig['one']
    
    def gaussian(x,mu,sigma,amplitude):
        return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )
    
    params_spe, params_cov_spe = scipy.optimize.curve_fit(gaussian,bins[:len(areaVals_notrig['one'])],spe_spectrum,p0=[0.4e8,1e7,400])
    
    return params_spe[0]

def main():
    m = get_gain("/bluearc/storage/SBC-17-data/","20170707_7")
    print(m)


if __name__ == "__main__":
    main()