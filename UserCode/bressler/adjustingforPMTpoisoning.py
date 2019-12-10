#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:41:40 2019

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
        if trace[i]>std and np.abs(t[i]-pk_loc)<2e-7:
            tPulse.append(t[i])
            pulse.append(trace[i])
    return [pulse,tPulse]

def SBC_pulse_integrator_bressler(trace,dt):
    baseline = np.mean(trace[0:100])
    baseline_std = np.std(trace[0:100])
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
        integral_t0 = tPMT[integral_to_index]
        tstep = tPMT[1]=tPMT[0]
        intwindow = 100 #ns
        integral_tf_index = intwindow/tstep + integral_t0_index
        p,t = get_pulse(trace,tPMT,integral_t0,baseline_std)
        ret = scipy.integrate.trapz(p,t)*CONVERSION_TO_CHARGE
        #ret=scipy.integral.trapz(trace[integral_t0_index:integral_tf_index],
                                 #tPMT[integral_t0_index:integral_tf_index])*CONVERSION_TO_CHARGE
    """
    if ret != None and Npeaks == 1 and a < 1e-7:   
        
        plt.figure()
        plt.plot(tPMT,trace)
        plt.plot(tPulse,pulse)
        plt.show()
    """     
    return [ret,Npeaks,totIntegral,pk_times]
    
def main():
    run = '20170706_1'
    runpath = "/bluearc/storage/SBC-17-data/"+run+'/'
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    allTraces = []
    totalAreas = []
    totareaofalltraces=[]

    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,0)
    tr = e["PMTtraces"]
    trac = tr["traces"]
    dt = tr["dt"]
    plotted = False
    for i in range(len(trac)):
        trace = np.fabs(trac[i][0])
        b = np.mean(trace[0:100])
        dt_tr = dt[i][0]
        tPMT = np.arange(len(trace))*dt_tr
        totareaofalltraces.append(total_area(trace-b,tPMT))
        integral_t0 = tPMT[np.argmax(np.diff(trace)>4)]
        if integral_t0 < 1e-7:
            if not plotted:
                print(integral_t0)
                plt.figure()
                plt.plot(tPMT,trace)
                plt.plot([min(tPMT),max(tPMT)],[b,b])
                plt.show
                plt.figure()
                plt.plot(tPMT[:-1],np.diff(trace))
                plt.show
                plotted = True

if __name__ == "__main__":
    main()
