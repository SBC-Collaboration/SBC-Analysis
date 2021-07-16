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
from pulse_integrator import stitchTraces
from pulse_integrator import SBC_pulse_integrator_bressler
from pulse_integrator import total_area
from random import random
import gc
from gaincalc import get_gain
import canyoufilterpmttraces as sbcfilter

def spectrum(datapath, run, forcebins = False):
    """Returns the spectrum and histogram bins for the scintillation of all events in a run"""
    runpath = datapath + '/' + run
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    allTraces = []
    total_time = 0
    pulses = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    times = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    #camextratime = 25e-6
    for event in events:
        if int(event)> 3:
            break
        print(event)
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
        if e["slowDAQ"]["loaded"]:
            #print(e["fastDAQ"].keys())
            cgate = e["fastDAQ"]["CAMgate"]
            #dcam = np.diff(cgate)
            fdt = e['fastDAQ']['time']
            #camOffTimes = np.sort(np.array([fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]))
            
            #camOnTimes = np.sort(np.array([fdt[i] for i in range(len(dcam)) if dcam[i] < 0.5]))
            fddt = fdt[1]-fdt[0]
            tfast = fdt[-1]-fdt[0]
            LED_on = [fdt[i] for i in range(len(cgate)) if cgate[i]<-0.5]
            blockedFraction = ((len(LED_on)*fddt))/tfast
            print(blockedFraction)
            tr = e["PMTtraces"]
            trac = tr["traces"]
            dt = tr["dt"]
            #event_time = (tr['t0_sec'][-1]+tr['t0_frac'][-1]-tr['t0_sec'][0] - tr['t0_frac'][0])[0]
            event_time = (((e["slowDAQ"]["elapsed_time"][-1]-e["slowDAQ"]["elapsed_time"][0]))*(1-blockedFraction))
            #print(event_time)
            total_time += event_time

            #f,axes = plt.subplots(1,5,sharey=True)
            #f.suptitle(runpath+"/"+str(event))
            #pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
            #d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
            #pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
            #tracetimes = pmttracetime - pmtalign
            #camoffindex = 0
            #camonindex = 0
            for i in range(len(trac)):
                #print(i)
                """
                thistracetime = tracetimes[i]
                
                #nearestcamoff = min(camOffTimes, key=lambda x:abs(x-thistracetime))
                #nearestcamon = min(camOnTimes, key=lambda x:abs(x-thistracetime))
                print(camOffTimes[camoffindex])
                print(thistracetime)
                if thistracetime > camOffTimes[camoffindex]:
                    camoffindex += 1
                if thistracetime > camOnTimes[camonindex]:
                     camonindex += 1 
                if camoffindex<len(camOffTimes)-1:
                    if abs(camOffTimes[camoffindex]-thistracetime)<camextratime:
                        print('excluding a trace near a camera off')
                        continue
                if camonindex<len(camOnTimes)-1:
                    if abs(camOnTimes[camonindex]-thistracetime)<camextratime:
                        print('excluding a trace near a camera on')
                        continue
                """
                trace = np.fabs(trac[i][0])
                if max(trace) == 128:
                    trace = stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][i][1]))
                dt_tr = dt[i][0]

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
                    allTraces.append(a)
                """
                #if a != None:
                if isZero:
                    if j < 5:
                        if isNegative:
                            if random() >0:
                                print(runpath+"/"+str(event)+" pmt trace "+str(i))
                                tPMT = np.arange(len(trace))*dt_tr
                                axes[j].plot(tPMT,trace,lw=3)
                                axes[j].set_xlabel("time (s)",fontsize=25)
                                axes[j].set_ylabel("PMT response (ADC)",fontsize=25)
                                j+=1
               
    
                plt.show
                """
                pulses[number].append(a)
        gc.collect()
            
            
    for k in pulses:
        pulses[k] = [x for x in pulses[k] if x != None]
    
    allTraces = [x for x in allTraces if x != None]
    
    plt.figure()

    Nbins = int(np.floor(np.sqrt(len(allTraces))))
    allvals, bins, _ = plt.hist(allTraces,Nbins,label='all traces')
    
    areaVals = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    for k in pulses:
        if k != 'other':
            areaVals[k], _, _ = plt.hist(pulses[k],bins,histtype = 'step',
                    linewidth = 3,label= k+' hits')
    plt.legend(fontsize=12)
    plt.show()    
    #spe_spectrum = areaVals['one']
    
    #def gaussian(x,mu,sigma,amplitude):
    #    return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )
    
    #params_spe, params_cov_spe = scipy.optimize.curve_fit(gaussian,bins[:len(areaVals['one'])],
    #                                                                    spe_spectrum,
    #                                                                    p0=[0.4e8,1e7,40])
    #params_twohits, params_cov_twohits = scipy.optimize.curve_fit(gaussian,
    #                                                              bins[:len(areaVals['two'])],
    #                                                              areaVals['two'],
    #                                                              p0=[0.8e8,1e7,10])
    #mu_spe = params_spe[0]
    #mu_2 = params_twohits[0] - mu_spe
    #print(mu_spe)
    #print(mu_2)
    
    #mu_avg = (mu_spe + mu_2)*0.5
    #mu_avg = get_gain(datapath,run)
    mu_avg = 1e7
    print(mu_avg)

    
    plt.figure()
    plt.grid(True)
    if isinstance(forcebins,np.ndarray):
        bins=forcebins
        fullspect,_,_=plt.hist([t/mu_avg for t in allTraces],
                                  forcebins,label='all traces')
        
    else:
        fullspect,bins,_=plt.hist([t/mu_avg for t in allTraces],
                                  int(np.floor(np.sqrt(len(allTraces)))),label='all traces')
        
    #print(bins)
    plt.yscale('log')
    plt.xlabel('phe based on a gain of '+str(mu_avg)+' electrons per phe')
    plt.legend()
    plt.show
    print(sum(fullspect)/total_time)
    print("The Total Exposure Time of run "+str(runpath)+ " was "+str(total_time)+" Seconds")
    print("The overall PMT trigger rate was " + str(len(allTraces)/total_time)+ "Hz")
    return [fullspect,bins,total_time]

def main():
    bgrun = '20170704_0' 
    bgrunpath = "/bluearc/storage/SBC-17-data/"
    bins = np.arange(-1,200)
    #bins = np.array([-1,0,0.2,0.4,0.6,0.8,1,1.2,1.3,1.4,1.6,1.8,2,3])
    bgspect,_,bgtime = spectrum(bgrunpath,bgrun,forcebins=bins)
    bgspect /= bgtime
    ibg = scipy.integrate.trapz(bgspect)
    
    print("trigger rate: "+str(ibg)+" Hz")
    bgspect /= bgtime
    plt.figure()
    plt.plot(bins[:len(bgspect)],bgspect,'g',linewidth=4,label='bi-207 sept 24')
    plt.xlabel('Photoelectrons',fontsize=25)
    plt.ylabel('Rate [Hz/phe]',fontsize=25)
    plt.yscale('log')
    plt.legend(fontsize=18)
    plt.show
    

if __name__ == "__main__":
    main()
