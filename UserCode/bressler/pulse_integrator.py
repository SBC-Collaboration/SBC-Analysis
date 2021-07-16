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
import random

# Global variable because I'm a physicist not a developer:
CONVERSION_TO_CHARGE = (125.0/128)*(1/50.0)*(1/1000.0)*(1/(1.602e-19))


def total_area(trace,t):
    """ Gets total area of trace with time array t"""
    return scipy.integrate.trapz(trace,x=t)*CONVERSION_TO_CHARGE

def get_pulse(trace, t, dt, locale, pk_loc, std):
    tPulse = []
    pulse = []
    tracemaxi = list(trace).index(max(trace))
    #print("tracemaxi: %d"%tracemaxi)
    for i in range(len(t)):
        if trace[i] < std and i > tracemaxi:
            break
        #print(np.abs(t[i]-pk_loc))
        elif trace[i]>=std and np.fabs(t[i]-pk_loc) <= locale:
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
    trace = trace[0:-100]

    pk_ind = scipy.signal.find_peaks(trace,5)
    #print(pk_ind)
    pk_times = [pk*dt for pk in pk_ind[0]]
    pk_vals = [trace[k] for k in pk_ind[0]]
    Npeaks = len(pk_vals)
    tPMT = np.arange(len(trace))*dt
    

    totIntegral = total_area(trace,tPMT)
    if Npeaks == 1:
        [pulse,tPulse] = get_pulse(trace, tPMT, dt, 0.5e-7, pk_times[0], baseline_std)
        ret = 0
        startind = 0
        
        for j in range(len(tPulse)-1):
            dist = tPulse[j+1] - tPulse[j]
            if dist>dt+1e-9:
                #print("break in t array at %d"%j)
                ret += scipy.integrate.trapz(pulse[startind:j+1],tPulse[startind:j+1])*CONVERSION_TO_CHARGE
                #print("ret inside pulse_integrator: %f"%ret)
                
                startind = j+1
            elif j == len(tPulse) - 2:
                #print("end of pulse condition, j = %d, t = %e"%(j,tPulse[j]))
                ret += scipy.integrate.trapz(pulse[startind:j+1],tPulse[startind:j+1])*CONVERSION_TO_CHARGE
    
        
    elif Npeaks ==0:
        [pulse,tPulse] = get_pulse(trace, tPMT, dt, 0.5e-7, 200*dt, baseline_std)
        ret = 0
        startind = 0
        #print(t)
        for j in range(len(tPulse)-1):
            dist = tPulse[j+1] - tPulse[j]
            if dist>dt+1e-9:
                #print("break in t array at %d"%j)
                ret += scipy.integrate.trapz(pulse[startind:j+1],tPulse[startind:j+1])*CONVERSION_TO_CHARGE
                #print("ret inside pulse_integrator: %f"%ret)
                startind = j+1
            elif j == len(tPulse) - 2:
                #print(j)
                ret += scipy.integrate.trapz(pulse[startind:],tPulse[startind:])*CONVERSION_TO_CHARGE
        
    elif Npeaks == 2:
        if np.abs(pk_times[0]-pk_times[1])>=2e-7:
            [firstPulse, tFirstPulse] = get_pulse(trace, tPMT, dt, 0.5e-7, pk_times[0], baseline_std)
            ret = 0
            startind = 0
            #print(t)
            for j in range(len(tFirstPulse)-1):
                dist = tFirstPulse[j+1] - tFirstPulse[j]
                if dist>dt + 1e-9:
                    #print("break in t array at %d"%j)
                    ret += scipy.integrate.trapz(firstPulse[startind:j+1],tFirstPulse[startind:j+1])*CONVERSION_TO_CHARGE
                    #print("ret inside pulse_integrator: %f"%ret)
                    startind = j+1
                elif j == len(tFirstPulse) - 2:
                    #print(j)
                    ret += scipy.integrate.trapz(firstPulse[startind:],tFirstPulse[startind:])*CONVERSION_TO_CHARGE
            [secondPulse, tSecondPulse] = get_pulse(trace,tPMT,dt, 0.5e-7, pk_times[1],baseline_std)
            
            startind = 0
            #print(t)
            for j in range(len(tSecondPulse)-1):
                dist = tSecondPulse[j+1] - tSecondPulse[j]
                if dist>dt+1e-9:
                    #print("break in t array at %d"%j)
                    ret += scipy.integrate.trapz(secondPulse[startind:j+1],tSecondPulse[startind:j+1])*CONVERSION_TO_CHARGE
                    #print("ret inside pulse_integrator: %f"%ret)
                    startind = j+1
                elif j == len(tSecondPulse) - 2:
                    #print(j)
                    ret += scipy.integrate.trapz(secondPulse[startind:],tSecondPulse[startind:])*CONVERSION_TO_CHARGE
            """
            if randrange(100) == 1 :
                plt.figure()
                #plt.title("baseline=%s"%str(baseline_std))
                plt.plot(tPMT,trace)
                plt.plot(tFirstPulse,firstPulse,linewidth=3)
                plt.plot(tSecondPulse,secondPulse,linewidth=3)
                plt.show
            """
            
        else:
            #print('-1')
            Npeaks = -1
            integral_t0_index = np.argmax(np.diff(trace)>4)
            integral_t0 = tPMT[integral_t0_index]
            p,t = get_pulse(trace,tPMT,dt, 5e-7,integral_t0,baseline_std)
            ret = 0
            startind = 0
            #print(t)
            for j in range(len(t)-1):
                dist = t[j+1] - t[j]
                if dist>dt+1e-9:
                    #print("break in t array at %d"%j)
                    ret += scipy.integrate.trapz(p[startind:j+1],t[startind:j+1])*CONVERSION_TO_CHARGE
                    #print("ret inside pulse_integrator: %f"%ret)
                    startind = j+1
                elif j == len(t) - 2:
                    #print(j)
                    ret += scipy.integrate.trapz(p[startind:],t[startind:])*CONVERSION_TO_CHARGE
  

    else:

        integral_t0_index = np.argmax(np.diff(trace)>4)
        integral_t0 = tPMT[integral_t0_index]
        p,t = get_pulse(trace, tPMT, dt, 5e-7, integral_t0, baseline_std)
        ret = 0
        startind = 0
        #print(t)
        for j in range(len(t)-1):
            dist = t[j+1] - t[j]
            if dist>dt+1e-9:
                #print("break in t array at %d"%j)
                ret += scipy.integrate.trapz(p[startind:j+1],t[startind:j+1])*CONVERSION_TO_CHARGE
                #print("ret inside pulse_integrator: %f"%ret)
                startind = j+1
            elif j == len(t) - 2:
                #print(j)
                ret += scipy.integrate.trapz(p[startind:],t[startind:])*CONVERSION_TO_CHARGE

        
        """
        if random.random()<0.001:   
            
            plt.figure()
            plt.plot(tPMT,trace)
            plt.plot(t,p)
            plt.xlabel('time (s)')
            plt.ylabel('signal (ADC units)')
            plt.show()
       """
    return [ret,Npeaks,totIntegral,pk_times]
    
def main():
    run = '20170709_8'
    runpath = "/bluearc/storage/SBC-17-data/"+run+'/'
    event = 0
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
    tr = e["PMTtraces"]
    trac = tr["traces"]
    dt = tr["dt"]
    for i in range(len(trac)):
        trace = np.fabs(trac[i][0])
        rawtrace = trac[i][0]
        stitched = False
        if max(trace) == 128:
            trace = stitchTraces(trace, np.fabs(trac[i][1]))
            stitchedtrace = stitchTraces(np.fabs(trac[i][0]), np.fabs(trac[i][1]))
            stitched=True
        else:
            abstrace = np.fabs(trac[i][0])
        b = np.mean(trace[0:50])
        trace -= b
        if stitched:
            bgsubtracestitched = stitchedtrace - b
        else:
            bgsubtrace = abstrace - b

        bstd = np.std(trace[:50])
        trace = trace[:-100]
        

        dt_tr = dt[i][0]
        tPMT = dt_tr*range(len(trace))
        
        [a_desired, Npeaks, a_tot, pk_times] = SBC_pulse_integrator_bressler(trace, dt_tr)
        if Npeaks == 1:
            #print("trace %d is a single peak, area %f"%(i, a_desired))
            if random.random() < 1./500:
                pk_ind = scipy.signal.find_peaks(trace,5)

                pk_times = [pk*dt_tr for pk in pk_ind[0]]
                [pulse,tPulse] = get_pulse(trace, tPMT, dt_tr, 0.5e-7, pk_times[0], bstd)
                
                plt.figure()

                plt.subplot(2,2,1)
                plt.plot(dt_tr*range(len(rawtrace)), rawtrace, label = 'raw trace')
                plt.xlabel('time (s)')
                plt.ylabel('sig (adc)')
                plt.legend()
                
                plt.subplot(2,2,2)
                plt.plot(dt_tr*range(len(abstrace)), abstrace, label = 'trace after absolute value')
                plt.xlabel('time (s)')
                plt.legend()
                plt.ylabel('sig (adc)')
                
                plt.subplot(2,2,3)
                plt.plot(dt_tr*range(len(bgsubtrace)), bgsubtrace, label = 'trace after background subtraction')
                plt.xlabel('time (s)')
                plt.ylabel('sig (adc)')
                plt.legend()
                
                plt.subplot(2,2,4)
                plt.plot(dt_tr*range(len(trace)), trace, label = 'trace after all cleanup operations')
                startind=0
                for j in range(1,len(tPulse)):
                    dist = tPulse[j] - tPulse[j-1]
                    if dist>dt_tr+1e-9:
                        #print("break in t array at %d"%j)
                        plt.plot(tPulse[startind:j],pulse[startind:j], 'm', linewidth=7, alpha=0.7)
                        #print("ret inside pulse_integrator: %f"%ret)
                        
                        startind = j+1
                    elif j == len(tPulse) - 2:
                        #print("end of pulse condition, j = %d, t = %e"%(j,tPulse[j]))
                        plt.plot(tPulse[startind:j+1],pulse[startind:j+1], 'm', linewidth=7, alpha=0.7, label = 'integrated region')
                plt.plot(tPulse, pulse, 'g:', linewidth=4, label = 'pulse region, area of integrated parts %.2e'%a_desired)
                plt.scatter(pk_times[0], trace[pk_ind[0]], 50, 'b', label='peak')
                plt.legend()
                plt.xlabel('time (s)')
                plt.ylabel('sig (adc)')
                plt.show()
                
        elif Npeaks>2 and a_desired > 1e10 and random.random()<1./1:
            integral_t0_index = np.argmax(np.diff(trace)>4)
            integral_t0 = tPMT[integral_t0_index]
            p,t = get_pulse(trace, tPMT, dt_tr, 5e-7, integral_t0, bstd)
            
            plt.figure()

            plt.subplot(2,2,1)
            plt.plot(dt_tr*range(len(rawtrace)), rawtrace, label = 'raw trace, channel 0')
            plt.xlabel('time (s)')
            plt.ylabel('sig (adc)')
            plt.legend()
            
            plt.subplot(2,2,2)
            if stitched:
                plt.plot(dt_tr*range(len(stitchedtrace)), stitchedtrace, label = 'trace after absolute value and stitching')
            else:
                plt.plot(dt_tr*range(len(abstrace)), abstrace, label = 'trace after absolute value')
            plt.xlabel('time (s)')
            plt.legend()
            plt.ylabel('sig (adc)')
            
            plt.subplot(2,2,3)
            if stitched:
                plt.plot(dt_tr*range(len(bgsubtracestitched)), bgsubtracestitched, label = 'trace after background subtraction')
            else:
                plt.plot(dt_tr*range(len(bgsubtrace)), bgsubtrace, label = 'trace after background subtraction')
            plt.xlabel('time (s)')
            plt.ylabel('sig (adc)')
            plt.legend()
            
            plt.subplot(2,2,4)
            plt.plot(dt_tr*range(len(trace)), trace, label = 'trace after all cleanup operations')
            
            startind = 0
            #print(t)
            for j in range(len(t)-1):
                dist = t[j+1] - t[j]

                if dist>dt_tr+1e-9:
                    #print("break in t array at %d-%d, dist=%e"%(j,j+1,dist))
                    #print(t[startind:j+1])
                    #print(p[startind:j+1])
                    plt.plot(t[startind:j+1],p[startind:j+1], 'm', linewidth=7, alpha=0.7)
                    #print("ret inside pulse_integrator: %f"%ret)
                    startind = j+1
                    #print(startind)
                elif j == len(t) - 2:
                    #print(j)
                    plt.plot(t[startind:],p[startind:], 'm', linewidth=7, alpha=0.7, label='integrated region')
            plt.plot(t, p, 'g:', linewidth=4, label = 'pulse region, area of integrated parts %.2e'%a_desired)
            plt.plot([integral_t0, integral_t0], [min(trace), max(trace)], 'b', label=r'$t_0$')
            plt.legend()
            plt.xlabel('time (s)')
            plt.ylabel('sig (adc)')
            plt.show()
    """
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
   """ 

if __name__ == "__main__":
    main()
