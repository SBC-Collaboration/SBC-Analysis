#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 09:11:28 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
from random import randrange

# Global variable because I'm a physicist not a developer:
CONVERSION_TO_CHARGE = (125.0/128)*(1/50.0)*(1/1000.0)*(1/(1.602e-19))


def total_area(trace,t):
    """ Gets total area of trace with time array t"""
    return scipy.integrate.trapz(trace,x=t)*CONVERSION_TO_CHARGE

def get_pulse(trace,t,dt,pk_loc,std):
    tPulse = []
    pulse = []
    for i in range(len(t)):
        #print(np.abs(t[i]-pk_loc))
        if trace[i]>std and np.fabs(t[i]-pk_loc)<100e-9:
            tPulse.append(t[i])
            pulse.append(trace[i])
    return [pulse,tPulse]

def stitchTraces(ch0Trace,ch1Trace):
    j = list(ch0Trace).index(128)
    multiplier = 128/ch1Trace[j]
    ch1Trace = [x*multiplier for x in ch1Trace]

    for i in range(len(ch0Trace)):
        if ch0Trace[i] ==128:
            ch0Trace[i] = ch1Trace[i]
    return ch0Trace

def SBC_pulse_integrator_bressler(trace,dt):
    """
    takes:
        trace - flipped (and stitched, if desired) PMT trace
        dt    - time step
    returns: (as a list)
        ret         - area of pulse
        Npeaks      - number of peaks scipy found in the trace
        totIntegral - total area under trace
        pk_times    - times of the peaks scipy found
    """
    baseline = np.mean(trace[0:50])
    baseline_std = np.std(trace[0:50])
    trace = trace - baseline
    pk_ind = scipy.signal.find_peaks(trace,5)
    #print(pk_ind)
    pk_times = [pk*dt for pk in pk_ind[0]]
    pk_vals = [trace[k] for k in pk_ind[0]]
    Npeaks = len(pk_vals)
    tPMT = np.arange(len(trace))*dt

    totIntegral = total_area(trace,tPMT)
    if Npeaks == 1:
        [pulse,tPulse] = get_pulse(trace, tPMT, dt, pk_times[0],baseline_std)
        a = scipy.integrate.trapz(pulse,tPulse)*CONVERSION_TO_CHARGE
        ret = a
    elif Npeaks ==0:
        [pulse,tPulse] = get_pulse(trace,tPMT,dt,200*dt,0)
        a = scipy.integrate.trapz(pulse,tPulse)*CONVERSION_TO_CHARGE
        ret = a
    elif Npeaks == 2:
        if np.abs(pk_times[0]-pk_times[1])>=2e-7:
            [firstPulse, tFirstPulse] = get_pulse(trace,tPMT,dt,pk_times[0],baseline_std)
            [secondPulse, tSecondPulse] = get_pulse(trace,tPMT,dt,pk_times[1],baseline_std)
            a = scipy.integrate.trapz(firstPulse,tFirstPulse)*CONVERSION_TO_CHARGE + scipy.integrate.trapz(secondPulse,tSecondPulse)*CONVERSION_TO_CHARGE
            ret = a
        else:
            ret = None
    elif Npeaks == 3:
        if min([np.abs(pk_times[0] - pk_times[1]), np.abs(pk_times[0]-pk_times[2]),
                np.abs(pk_times[2]-pk_times[1])]) >= 2e-7:
            [firstPulse, tFirstPulse] = get_pulse(trace,tPMT,dt,pk_times[0],baseline_std)
            [secondPulse, tSecondPulse] = get_pulse(trace,tPMT,dt,pk_times[1],baseline_std)
            [thirdPulse, tThirdPulse] = get_pulse(trace,tPMT,dt,pk_times[2],baseline_std)
            a = scipy.integrate.trapz(firstPulse,tFirstPulse)*CONVERSION_TO_CHARGE + scipy.integrate.trapz(secondPulse,
                                     tSecondPulse)*CONVERSION_TO_CHARGE + scipy.integrate.trapz(thirdPulse,tThirdPulse)*CONVERSION_TO_CHARGE
            ret = a
        else:
            ret = None
            
    else:

        integral_t0_index = np.argmax(np.diff(trace)>4)
        integral_t0 = tPMT[integral_t0_index]
        p,t = get_pulse(trace,tPMT,dt,integral_t0,baseline_std)
        ret = scipy.integrate.trapz(p,t)*CONVERSION_TO_CHARGE
        """
        if randrange(10000) == 2:
            plt.figure()
            #plt.title("baseline=%s"%str(baseline_std))
            plt.plot(tPMT,trace)
            plt.plot(t,p,linewidth=3)
            plt.plot([integral_t0, integral_t0],[0,10],'r',linewidth=5)
            plt.show
        """
        
    """
    if ret != None and Npeaks == 1 and a < 1e-7:   
        
        plt.figure()
        plt.plot(tPMT,trace)
        plt.plot(tPulse,pulse)
        plt.show()
    """     
    return [ret,Npeaks,totIntegral,pk_times]
    
def main():
    run = '20170706_4'
    runpath = "/bluearc/storage/SBC-17-data/"+run+'/'
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    allTraces = []
    totalAreas = []
    totareaofalltraces=[]
    selftrig_pulses = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    notrig_pulses = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    for event in events:
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
        tr = e["PMTtraces"]
        trac = tr["traces"]
        dt = tr["dt"]
        for i in range(len(trac)):
            trace = np.fabs(trac[i][0])
            b = np.mean(trace[0:100])
            NIMtrace = trac[i][1]
            # Determine whether the NIM module would've seen a trigger:
            if min(NIMtrace) < -30:
                NTrig = True
            else:
                NTrig = False
            # get the time step, assuming it's constant
            dt_tr = dt[i][0]
            tPMT = np.arange(len(trace))*dt_tr
            totareaofalltraces.append(total_area(trace-b,tPMT))

            # populate dictionaries arrays based on how many pulses there were
            if NTrig:
                [a_trig,n_trig,totInt,pktimes] = SBC_pulse_integrator_bressler(trace,dt_tr)
                if n_trig == 0:
                    number = 'zero'
                elif n_trig == 1:
                    number = 'one'
                elif n_trig == 2:
                    number = 'two'
                elif n_trig == 3:
                    number = 'three'
                else:
                    number = 'other'
                selftrig_pulses[number].append(a_trig)
                allTraces.append(a_trig)

            else:
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
                
            totalAreas.append(totInt)
    for k in selftrig_pulses:
        selftrig_pulses[k] = [x for x in selftrig_pulses[k] if x != None]
    for k in notrig_pulses:
        notrig_pulses[k] = [x for x in notrig_pulses[k] if x != None]
    
    allTraces = [x for x in allTraces if x != None]
    totalAreas = [x for x in totalAreas if x != None]
    
    plt.figure()
    plt.grid(True)
    plt.hist([t/4e7 for t in totareaofalltraces],int(np.floor(np.sqrt(len(totareaofalltraces)))))
    plt.yscale('log')
    plt.xlabel('phe based on a gain of 4e7 electrons per phe')
    plt.show
    
    plt.figure()
    Nbins = int(np.floor(np.sqrt(len(allTraces))))
    allvals, bins, _ = plt.hist(allTraces,Nbins,label='all traces')
    
    areaVals_trig = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    areaVals_notrig = {'zero': [], 'one': [], 'two': [], 'three': [], 'other': []}
    
    for k in selftrig_pulses:
        areaVals_trig[k], _, _ = plt.hist(selftrig_pulses[k],bins,histtype = 'step',linewidth = 3,label='trig '+k+' hits')
        areaVals_notrig[k], _, _ = plt.hist(notrig_pulses[k],bins,histtype = 'step',linewidth = 3,label='no trig '+k+' hits')
        
    spe_spectrum = areaVals_trig['one']+areaVals_notrig['one']
    
    def gaussian(x,mu,sigma,amplitude):
        return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )
    
    params_spe, params_cov_spe = scipy.optimize.curve_fit(gaussian,bins[:len(areaVals_trig['one'])],spe_spectrum,p0=[0.4e8,1e7,400])
    params_twohits, params_cov_twohits = scipy.optimize.curve_fit(gaussian,bins[:len(areaVals_trig['two'])],areaVals_trig['two'],p0=[0.8e8,1e7,10])
    print('spe fit:')
    print('mu = '+str(params_spe[0]/1e7)+'*10^7')
    print('sigma = '+str(params_spe[1]/1e7)+'*10^7')
    print('\n')
    print('two-hit fit:')
    print('mu = '+str(params_twohits[0]/1e7)+'*10^7')
    print('sigma = '+str(params_twohits[1]/1e7)+'*10^7')
    plt.plot(bins[:len(areaVals_trig['two'])],gaussian(bins[:len(areaVals_trig['two'])],params_twohits[0],params_twohits[1],params_twohits[2]),
                  color='r',linewidth=5,label='Fit to two peak distribution, mu='+str(params_twohits[0]))
    plt.plot(bins[:len(areaVals_trig['one'])],gaussian(bins[:len(areaVals_trig['one'])],params_spe[0],params_spe[1],params_spe[2]),
                  color='m',linewidth=5,label='Fit to spe  distribution, mu='+str(params_spe[0]))

    plt.legend()
    plt.title(run)
    plt.xlabel('Charge (electrons)')
    plt.yscale('log')
    plt.ylim([0.5, 1e3])
    plt.rcParams.update({'font.size':18})
    plt.show
    
    plt.figure()
    plt.plot(bins[:len(areaVals_trig['one'])],spe_spectrum)
    plt.plot(bins[:len(areaVals_trig['one'])],gaussian(bins[:len(areaVals_trig['one'])],params_spe[0],params_spe[1],params_spe[2]),
                  color='r',linewidth=5)
    plt.title('spe spectrum with gaussian fit')
    plt.show
    
    
    effXvals = bins[:len(areaVals_trig['one'])]
    effXvals_NIM = bins[:len(areaVals_notrig['one'])]
    effic_byNpeak = np.divide(spe_spectrum,allvals)
    effic_byNIM = np.divide(areaVals_trig['one'],allvals)
    effic_byNpeak[np.isnan(effic_byNpeak)]=1#float('+inf')
    effic_byNIM[np.isnan(effic_byNIM)]=1#float('+inf')
    #effXvals_NIM = effXvals_NIM[effic_byNIM<float('+inf')]
    #effic_byNIM=effic_byNIM[effic_byNIM<float('+inf')]
    
    def functn(x,a,b):
        return scipy.stats.norm.cdf(x,a,b)
    gausspec = gaussian(bins[:len(areaVals_trig['one'])],params_spe[0],params_spe[1],params_spe[2])
    eff_params, eff_params_cov = scipy.optimize.curve_fit(functn,bins[:len(effic_byNIM)],effic_byNIM)
    #eff_dens = np.multiply(effic_byNIM,spe_spectrum)
    eff_dens = np.multiply(effic_byNIM,gausspec)
    numerator = scipy.integrate.trapz(eff_dens,effXvals)
    #denominator = scipy.integrate.trapz(spe_spectrum,effXvals)
    denominator = scipy.integrate.trapz(gaussian(bins[:len(areaVals_trig['one'])],params_spe[0],params_spe[1],params_spe[2]),bins[:len(areaVals_trig['one'])])
    print('efficiency = '+str(numerator/denominator))
    
    plt.figure()
    #plt.plot(bins[:len(effic_byNIM)],functn(bins[:len(effic_byNIM)],params[0],params[1]),color='r')
    #plt.text(40,.75,"mu = "+str(params[0]),fontsize=15)
    #plt.text(40,.5,"sigma = "+str(params[1]),fontsize=15)
    plt.plot(effXvals,effic_byNpeak,linewidth=3)
    plt.plot(effXvals_NIM,effic_byNIM,linewidth=3)
    plt.xlabel('area',fontsize=18)
    plt.ylabel('efficiency',fontsize=18)
    plt.show
    

if __name__ == "__main__":
    main()
