#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 10:24:55 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pulse_integrator as pi
from gaincalc import get_gain


def NIM_efficiency_and_plot(V, VwithNIM, A, AwithNIM, Anogain, AnogainNIM, title_string, smallPulses, smallPulseBins):
    plt.figure()
    vvals, vbins, _= plt.hist(np.asarray(V),126,color='r',histtype = 'step',label='All Traces, N=%d'%len(V))
    vnimvals, _, _ = plt.hist(np.asarray(VwithNIM),bins=vbins,color='b',histtype='step',label='Traces with NIM signal, N=%d'%len(VwithNIM))
    plt.xlabel('signal max [ADC]', fontsize=18)
    plt.legend(fontsize=15)
    plt.title(title_string)
    plt.ylabel('Count', fontsize=18)
    plt.grid()
    plt.yscale('log')
    plt.show
    
    plt.figure()
    Anogainvals, Abins, _ = plt.hist(Anogain, 224, histtype='step', color='k', linewidth=2, label='all data')

    Abincenters = [0.5*(Abins[i+1]+Abins[i]) for i in range(len(Abins)-1)]
        
    def gaussian(x, mu, sigma, amplitude):
        return (amplitude/(np.sqrt(2*np.pi)*sigma))*np.exp(-0.5*((x-mu)/sigma)**2)
    
    onepestart = 0.5e7
    onepeend = 2e7
    afit = [Anogainvals[i] for i in range(len(Anogainvals)) if Abincenters[i] > onepestart and Abincenters[i] < onepeend]
    bincfit = [x for x in Abincenters if x > onepestart and x < onepeend]
    startpoint = [1e7,2e7,2000]
    plt.scatter(bincfit, afit, 15, color='b', label='fitted data', zorder=2.5)

    p, pcov = scipy.optimize.curve_fit(gaussian,bincfit,
                                       afit,p0=startpoint)
    print(p)
    x = np.arange(1,5e7,0.01e7)
    #x=bincfit
    plt.plot(x,[gaussian(y,p[0],p[1],p[2]) for y in x], 'b:', linewidth=1, label=r'fit, $\mu=$%.2e'%p[0])
    plt.hist(smallPulses, smallPulseBins, linewidth=2, histtype='step', color='r', label='single peak pulses')
    plt.xlim([0,10e7])
    plt.legend(fontsize=15)
    plt.xlabel('Pulse area', fontsize=18)
    plt.ylabel('counts', fontsize=18)
    plt.grid()
    plt.yscale('log')
    plt.title(title_string)
    plt.show()
    
    plt.figure()
    avals, bins, _= plt.hist(A,244,histtype = 'step',label='All Traces, N=%d'%len(A))
    animvals, _, _ = plt.hist(AwithNIM,bins=bins,histtype='step',label='Traces with NIM signal, N=%d'%len(AwithNIM))
    plt.xlabel('signal [phe]')
    plt.legend()
    plt.title(title_string)
    plt.show
    
    plt.figure()
    AnogainNIMvals, _, _ = plt.hist(AnogainNIM, Abins)
    plt.show()

    Aadjustedgain = [x/p[0] for x in Anogain]
    ANIMadjustedgain = [x/p[0] for x in AnogainNIM]
    plt.figure()
    plt.hist(A,bins=bins,histtype='step',label='All Traces, "normal" gain calculation')
    aadjustedvals, _, _ = plt.hist(Aadjustedgain, bins=bins, histtype='step', label='All Traces, gain calculated from all traces fit')
    plt.hist(AwithNIM, bins=bins, histtype='step', label='traces with NIM signal, old gain')
    animadjustedvals, _, _ = plt.hist(ANIMadjustedgain, bins=bins, histtype='step', label='traces with NIM, new gain')
    plt.xlabel('phe')
    plt.ylabel('count')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.title(title_string)
    plt.show()
    
    plt.figure()
    plt.hist(Aadjustedgain, bins=bins, histtype='step', label='All Traces, gain calculated from all traces fit, N=%d'%len(Aadjustedgain))
    plt.hist(ANIMadjustedgain, bins=bins, histtype='step', label='traces with NIM, new gain, N=%d'%len(ANIMadjustedgain))
    plt.xlabel('phe')
    plt.ylabel('count')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.title(title_string)
    plt.show()
    bincenters = [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1)]

    
    animadjustedvals = animadjustedvals[aadjustedvals>0]
    aadjustedvals = aadjustedvals[aadjustedvals>0]
    afrac = np.divide(animadjustedvals, aadjustedvals)
    afrac[np.isnan(afrac)] = float("+inf")
    afrac = afrac[afrac<float("+inf")]

    AnogainNIMvals = AnogainNIMvals[Anogainvals>0]
    Anogainvals = Anogainvals[Anogainvals>0]
    anogainfrac = np.divide(AnogainNIMvals, Anogainvals)
    anogainfrac[np.isnan(anogainfrac)] = float("+inf")
    anogainfrac = anogainfrac[anogainfrac<float("+inf")]

    vnimvals = vnimvals[vvals>0]
    vvals = vvals[vvals>0]
    
    perc = np.divide(vnimvals,vvals)
    perc[np.isnan(perc)]=float('+inf')
    perc=perc[perc<float('+inf')]
    
    def functn(x,a,b):
        return scipy.stats.norm.cdf(x,a,b)
    
    params, params_cov = scipy.optimize.curve_fit(functn,vbins[:len(perc)],perc,p0=[1,1])
    params_pheefficiency, params_cov_pheefficiency = scipy.optimize.curve_fit(functn, 
                                                                              bincenters[:len(afrac)], afrac, p0=[1,2])
    params_areaefficiency, params_cov_areaefficiency = scipy.optimize.curve_fit(functn, 
                                                                              Abincenters[:len(anogainfrac)], anogainfrac, p0=[0.2e7,2e6])
    
    AwithNIMandEfficiency = [animadjustedvals[i]/functn(bincenters[i], 
                             params_pheefficiency[0], params_pheefficiency[1]) for i in range(len(afrac))]
    AnogainNIMandEfficiency = [AnogainNIMvals[i]/functn(Abincenters[i], 
                             params_areaefficiency[0], params_areaefficiency[1]) for i in range(len(anogainfrac))]
    
    plt.figure()
    plt.hist(Aadjustedgain, bins=bins, histtype='step',label='all traces, adjusted gain')
    plt.hist(ANIMadjustedgain, bins=bins, histtype='step', label='traces with NIM, adjusted gain')
    plt.scatter(bincenters[:len(afrac)],AwithNIMandEfficiency, label='traces with NIM, divided by efficiency')
    plt.legend()
    plt.xlabel('phe')
    plt.ylabel('count')
    plt.yscale('log')
    plt.title(title_string)
    plt.grid()
    plt.show()
    
    plt.figure()
    plt.hist(Anogain, bins=Abins, histtype='step',label='all traces, no gain')
    plt.hist(AnogainNIM, bins=Abins, histtype='step', label='traces with NIM, no gain')
    plt.scatter(Abincenters[:len(anogainfrac)],AnogainNIMandEfficiency, label='traces with NIM, divided by efficiency')
    #plt.scatter(Abincenters[:len(Anogainvals)], Anogainvals)
    #plt.scatter(Abincenters[:len(AnogainNIMvals)], AnogainNIMvals)
    print(params_areaefficiency)
    plt.legend()
    plt.xlabel('electrons')
    plt.ylabel('count')
    plt.yscale('log')
    plt.title(title_string)
    plt.grid()
    plt.show()
    
    plt.figure()
    plt.scatter(bincenters[:len(afrac)],afrac)
    plt.plot(bincenters[:len(afrac)], functn(bincenters[:len(afrac)], params_pheefficiency[0], params_pheefficiency[1]), color='r')
    plt.xlabel('phe')
    plt.ylabel('efficiency')
    plt.text(2.5,.75,"mu = %.2f"%(params_pheefficiency[0]),fontsize=15)
    plt.text(2.5,.5,"sigma = %.2f"%(params_pheefficiency[1]),fontsize=15)
    plt.show()
    
    plt.figure()
    plt.scatter(vbins[:len(perc)],perc)
    plt.plot(vbins[:len(perc)],functn(vbins[:len(perc)],params[0],params[1]),color='r')
    plt.text(40,.75,"mu = %.2f"%(params[0]),fontsize=15)
    plt.text(40,.5,"sigma = %.2f"%(params[1]),fontsize=15)
    plt.xlabel('signal max [ADC]')
    plt.ylabel('efficiency')
    plt.title(title_string)
    plt.show()

    
def main():
    run = '20170706_4'
    datapath = "/bluearc/storage/SBC-17-data/"
    runpath = datapath + run
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    [m, smallpulses, smallpulsebins] = get_gain(datapath,run, sendPulses=True)
    print("gain calculated with all small pulses: %.2e"%m)
    V = []
    VwithNIM = []
    A = []
    AwithNIM = []
    Anogain = []
    AnogainNIM = []
    for event in events:
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
        tr = e["PMTtraces"]
        trac = tr["traces"]
        dt = tr["dt"]
        for i in range(len(trac)):
            trace = np.fabs(trac[i][0])
            nimpulse = trac[i][1]
            b = np.mean(trace[0:50])
            trace -= b
            dt_tr = dt[i][0]
            V.append(max(trace))
            A.append(pi.SBC_pulse_integrator_bressler(trace,dt_tr)[0]/m)
            Anogain.append(pi.SBC_pulse_integrator_bressler(trace,dt_tr)[0])
    
            if min(nimpulse) < -30:
                VwithNIM.append(max(trace))
                AwithNIM.append(pi.SBC_pulse_integrator_bressler(trace,dt_tr)[0]/m)
                AnogainNIM.append(pi.SBC_pulse_integrator_bressler(trace,dt_tr)[0])
    NIM_efficiency_and_plot(V, VwithNIM, A, AwithNIM, Anogain, AnogainNIM, '20170706_4', smallpulses, smallpulsebins)

    
if __name__=="__main__":
    main()
    

