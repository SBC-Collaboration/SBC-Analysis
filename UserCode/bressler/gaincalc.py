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
from pulse_integrator import SBC_pulse_integrator_bressler, get_pulse
from pulse_integrator import stitchTraces
import canyoufilterpmttraces as sbcfilter
import random

def get_gain(datapath,run,sendPulses = False):
    
    runpath = datapath+"/"+run+'/'
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    allTraces = []
    pulses = {'zero': [], 'one': [], 'other': []}
    optm = []
    aptm = []
    totInts = []
    for event in events:
        #plt.figure()
        if int(event) < 10:
            e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
    
            tr = e["PMTtraces"]
            trac = tr["traces"]
            dt = tr["dt"]
            for i in range(len(trac)):
                stitched = False
                trace = np.fabs(trac[i][0])
                # get the time step, assuming it's constant
                dt_tr = dt[i][0]
                
                if max(trace) == 128:
                    stitched = True
                    trace = stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][i][1]))
                #trace = sbcfilter.filterXeBCTrace(trace,dt_tr)
                b = np.mean(trace[0:50])
                trace -= b
                
    
                # populate dictionaries arrays based on how many pulses there were
                
                
                [a,n,totInt,pktimes] = SBC_pulse_integrator_bressler(trace,dt_tr)
                
                if n == 0:
                    number = 'zero'
                elif n == 1:
                    number = 'one'
                    #print('single pulse event %d trace %i area %.2f'%(int(event), i, a))
                    if stitched:
                        print('stitched and single pulse!')
                    optm.append(max(trace))
                    totInts.append(totInt)
                    trace = trace[:-100]
                    """
                    if random.random()<0.005:
                        pk_ind = scipy.signal.find_peaks(trace,5)
                        pk_times = [pk*dt_tr for pk in pk_ind[0]]
                        pk_vals = [trace[k] for k in pk_ind[0]]
                        tPMT = np.arange(len(trace))*dt_tr
    
                        p, tp = get_pulse(trace,tPMT,dt_tr, 0.5e-7,pk_times[0], np.std(trace[0:50]))
                        plt.figure()
                        plt.plot(dt_tr * range(len(trace)), trace, 'r')
                        
                        startind = 0
                        for j in range(1,len(tp)):
                            dist = tp[j] - tp[j-1]
                            if dist>dt_tr+1e-9:
                                plt.plot(tp[startind:j],p[startind:j],'k', linewidth=3)
                                startind = j
                            elif j == len(tp) - 2:
                                plt.plot(tp[startind:j+1], p[startind:j+1], 'k', linewidth=3)
                        plt.scatter(pk_times, pk_vals, 50)
                        #plt.title('low area trace')
                        plt.show()
                    
                    elif a>1e7 and random.random()<0.01:
                        pk_ind = scipy.signal.find_peaks(trace,5)
                        pk_times = [pk*dt_tr for pk in pk_ind[0]]
                        pk_vals = [trace[k] for k in pk_ind[0]]
                        tPMT = np.arange(len(trace))*dt_tr
    
                        p, tp = get_pulse(trace,tPMT,dt_tr, 0.5e-7,pk_times[0], np.std(trace[0:50]))
                        plt.figure()
                        plt.plot(dt_tr * range(len(trace)), trace, 'b')
                        startind = 0
                        for j in range(1,len(tp)):
                            dist = tp[j] - tp[j-1]
                            if dist>dt_tr:
                                plt.plot(tp[startind:j],p[startind:j],'k', linewidth=3)
                                
                                startind = j
                            elif j == len(tp) - 2:
                                plt.plot(tp[startind:j+1], p[startind:j+1], 'k', linewidth=3)
                        plt.scatter(pk_times, pk_vals, 50)
                        plt.title('high area trace')
                        plt.show()
                    """
                else:
                    number = "other"
    
                pulses[number].append(a)
                allTraces.append(a)
                aptm.append(max(trace))
        plt.show()
    for k in pulses:
        pulses[k] = [x for x in pulses[k] if x != None]
    
    allTraces = [x for x in allTraces if x != None]
    
    Nbins = int(np.floor(np.sqrt(len(allTraces))))
    allvals, bins, _ = plt.hist(allTraces,Nbins,label='all traces')
    areaVals = {'zero': [], 'one': [], 'other': []}
    
    Nbins = int(np.floor(np.sqrt(len(pulses["one"]))))
    areaVals["one"], bins = np.histogram(pulses["one"],Nbins)
        
    spe_spectrum = areaVals['one']
    
    bincenters = [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1)]
    spe_spectrum_efficiency = [spe_spectrum[i]/scipy.stats.norm.cdf(bincenters[i], 4820746.10225926, 1960279.07188277) for i in range(len(bincenters))]
    
    fitspect = [spe_spectrum[i] for i in range(len(spe_spectrum)) if bincenters[i]>0.2e7 and bincenters[i]<1e8]
    fitx = [x for x in bincenters if x>0.2e7 and x<1e8]
    
    def gaussian(x,mu,sigma,amplitude):
        return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )
    
    params_spe, params_cov_spe = scipy.optimize.curve_fit(gaussian,
                                                          fitx,
                                                          fitspect,p0=[1e7,4e7,3000])
    
    
    plt.figure()
    plt.hist(pulses["one"],bins,histtype='step',label='Single pulses')
    plt.scatter(bincenters,spe_spectrum_efficiency, label='after adjusting for efficiency')
    plt.plot(np.arange(0, 1e8, 1e6), 
             [gaussian(x, params_spe[0],params_spe[1],params_spe[2]) for x in np.arange(0, 1e8, 1e6)],
             label = r'Gaussian fit, $\mu=$%.2e, $\sigma=$%.2e'%(params_spe[0],params_spe[1]))
    plt.xlabel("Single Pulse Areas [electrons]")
    plt.legend(fontsize = 12)
    plt.title(run)
    plt.show

    
    plt.figure()
    h,x,y,p = plt.hist2d(optm,pulses["one"],100)
    plt.xlabel('trace max')
    plt.ylabel('area')
    plt.title('single pulses in '+run)
    plt.show()
    
    plt.figure()
    plt.hist(optm,100)
    plt.xlabel('trace max')
    plt.show()
    
    plt.figure()
    plt.hist(allTraces,100)
    plt.title("all traces area")
    plt.show()
    
    plt.figure()
    h,thesebins,_ = plt.hist(aptm, 100)
    plt.hist(optm,thesebins,histtype='step')
    plt.title('all traces max')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    if sendPulses:
        return [params_spe[0], pulses["one"], bins]
    else:
        return params_spe[0]
    
    
def main():
    m = get_gain("/bluearc/storage/SBC-17-data/",'20170719_0')
    print(m/1e7)


if __name__ == "__main__":
    main()