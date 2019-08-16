#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:28:31 2019

@author: bressler
"""


import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
from random import randrange
from pulse_integrator import SBC_pulse_integrator_bressler
from pulse_integrator import total_area

def main():
    run = '2017071_17'
    runpath = "/bluearc/storage/SBC-17-data/"+run+'/'
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    allTraces = []
    totalAreas = []
    totareaofalltraces=[]
    pulses = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    times = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
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
            tPMT = np.arange(len(trace))*dt_tr
            totareaofalltraces.append(total_area(trace-b,tPMT))

            # populate dictionaries arrays based on how many pulses there were
            [a,n,totInt,pktimes] = SBC_pulse_integrator_bressler(trace,dt_tr)
            if n == 0:
                number = 'zero'
                allTraces.append(a)
            elif n == 1:
                number = 'one'
                allTraces.append(a)
                times['one'].append(pktimes[0])
            elif n == 2:
                number = 'two'
                allTraces.append(a)
            elif n == 3:
                number = 'three'
                allTraces.append(a)
            else:
                number = 'other'
            pulses[number].append(a)
            
            
    for k in pulses:
        pulses[k] = [x for x in pulses[k] if x != None]
    
    allTraces = [x for x in allTraces if x != None]
    totalAreas = [x for x in totalAreas if x != None]
    
    plt.figure()
    Nbins = int(np.floor(np.sqrt(len(allTraces))))
    allvals, bins, _ = plt.hist(allTraces,Nbins,label='all traces')
    
    areaVals = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    
    for k in pulses:
        if k != 'other':
            areaVals[k], _, _ = plt.hist(pulses[k],bins,histtype = 'step',
                    linewidth = 3,label= k+' hits')
        
    spe_spectrum = areaVals['one']
    
    def gaussian(x,mu,sigma,amplitude):
        return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )
    
    params_spe, params_cov_spe = scipy.optimize.curve_fit(gaussian,bins[:len(areaVals['one'])],
                                                                        spe_spectrum,
                                                                        p0=[0.4e8,1e7,40])
    params_twohits, params_cov_twohits = scipy.optimize.curve_fit(gaussian,
                                                                  bins[:len(areaVals['two'])],
                                                                  areaVals['two'],
                                                                  p0=[0.8e8,1e7,10])
    mu_spe = params_spe[0]
    mu_2 = params_twohits[0]
    print('spe fit:')
    print('mu = '+str(mu_spe/1e7)+'*10^7')
    print('sigma = '+str(params_spe[1]/1e7)+'*10^7')
    print('\n')
    print('two-hit fit:')
    print('mu = '+str(mu_2/1e7)+'*10^7')
    print('sigma = '+str(params_twohits[1]/1e7)+'*10^7')
    plt.plot(bins[:len(areaVals['two'])],gaussian(bins[:len(areaVals['two'])],
             params_twohits[0],params_twohits[1],params_twohits[2]),
                  color='r',linewidth=5,
                  label='Fit to two peak distribution, mu='+str(params_twohits[0]))
    plt.plot(bins[:len(areaVals['one'])],gaussian(bins[:len(areaVals['one'])],params_spe[0],params_spe[1],params_spe[2]),
                  color='m',linewidth=5,
                  label='Fit to spe  distribution, mu='+str(params_spe[0]))

    plt.legend()
    plt.title(run)
    plt.xlabel('Charge (electrons)')
    plt.yscale('log')
    plt.rcParams.update({'font.size':18})
    plt.ylim([0.5,1e6])
    plt.show
    
    plt.figure()
    plt.plot(bins[:len(areaVals['one'])],spe_spectrum)
    plt.plot(bins[:len(areaVals['one'])],gaussian(bins[:len(areaVals['one'])],
             params_spe[0],params_spe[1],params_spe[2]),
                  color='r',linewidth=5)
    plt.title('spe spectrum with gaussian fit')
    plt.show
    
    plt.figure()
    plt.hist(times['one'],int(np.sqrt(len(times['one']))))
    plt.yscale('log')
    plt.show
    
    plt.figure()
    plt.grid(True)
    plt.hist([t/mu_spe for t in pulses['other']],int(np.floor(np.sqrt(len(pulses['other'])))))
    plt.yscale('log')
    plt.xlabel('phe based on a gain of '+str(mu_spe)+' electrons per phe')
    plt.show
    
    

if __name__ == "__main__":
    main()
