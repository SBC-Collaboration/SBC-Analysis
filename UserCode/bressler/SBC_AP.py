#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:27:33 2019

@author: bressler
"""

import SBCcode as sbc
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import runlistscatalogue as rlc
import gc

def XeBCAP(datadir,run,outdir,f_low,f_high):
    print("XeBCAP processing run %s"%run)
    runpath = datadir+'/'+run
    events = [evnt for evnt in os.listdir(runpath) if not os.path.isfile(os.path.join(runpath,evnt))]
    Nevents = len(events)
    #print(Nevents)
    outfilepath = outdir+'/'+'XeBCAP_'+run+'.txt'
    with open(outfilepath,'w') as fout:
        fout.write("event AP1subtracted AP1 AP1noise AP2subtracted AP2 AP2noise\n")
        for i in range(Nevents):
            #print(i)
            gc.collect()
            e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,i,'fastDAQ')
            if e["fastDAQ"]["loaded"]:
                fd = e['fastDAQ']
                t_full = fd['time']
                p1full = fd["Piezo1"]
                p2full = fd["Piezo2"]
                
                tstep = t_full[1]-t_full[0] #sample timestep, s
                Fs = 1/tstep #sampling frequency, Hz
                window_size = 0.02 #seconds
                
                runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
                acousticfilename = runreconpath+"AcousticAnalysis_%s.bin"%run
                a = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
                if i + 1 > len(a["bubble_t0"][:,0]):
                    print("i out of bounds, incomplete event")
                    continue
                at01 = a["bubble_t0"][i,0]
                at02 = a["bubble_t0"][i,1]
                
                p1_noise = [p1full[i] for i in range(len(p1full)) if t_full[i]-t_full[0] < window_size]
                p2_noise = [p2full[i] for i in range(len(p2full)) if t_full[i]-t_full[0] < window_size]
                
                p1 = [p1full[i] for i in range(len(p1full)) if t_full[i]>at01 and t_full[i]<at01+window_size]
                p2 = [p2full[i] for i in range(len(p2full)) if t_full[i]>at02 and t_full[i]<at02+window_size]
                N=len(p1)
                #print(len(p1_noise))
                #print(len(p1))
                
                if len(p1_noise) == len(p1):
                    #FFT stuff
                    fft1 = np.fft.rfft(p1)
                    fft1noise = np.fft.rfft(p1_noise)
            
                    fft2 = np.fft.rfft(p2)
                    fft2noise = np.fft.rfft(p2_noise)
                    
                    #Power Spectrum stuff
                    ps1 = (1/(N*Fs)) * (np.absolute(fft1)**2)
                    ps1noise = (1/(N*Fs)) * (np.absolute(fft1noise)**2)
                    ps2 = (1/(N*Fs)) * (np.absolute(fft2)**2)
                    ps2noise = (1/(N*Fs)) * (np.absolute(fft2noise)**2)        
                    freq = np.arange(0,Fs/2+Fs/N,Fs/N)
                    """
                    print("|ps1|="+str(len(ps1)))
                    print("|ps1noise|="+str(len(ps1noise)))
                    print("|ps2|="+str(len(ps2)))
                    print("|ps2noise|="+str(len(ps2noise)))
                    print("|freq|="+str(len(freq)))
                    """
        
                    #cut to the right frequencies
                    #print(len(freq))
                    #print(len(ps1))
                    #print(len(ps1noise))
                    ps1 = [ps1[i] for i in range(len(ps1)) if freq[i] > f_low and freq[i] < f_high]
                    ps1noise = [ps1noise[i] for i in range(len(ps1noise)) if freq[i] > f_low and freq[i] < f_high]
                    ps2 = [ps2[i] for i in range(len(ps2)) if freq[i] > f_low and freq[i] < f_high]
                    ps2noise = [ps2noise[i] for i in range(len(ps2noise)) if freq[i] > f_low and freq[i] < f_high]
                    freq = [x for x in freq if x > f_low and x < f_high]
            
                    #integrate
                    AP1 = scipy.integrate.trapz(ps1,freq)
                    AP1noise = scipy.integrate.trapz(ps1noise,freq)
                    
                    AP2 = scipy.integrate.trapz(ps2,freq)
                    AP2noise = scipy.integrate.trapz(ps2noise,freq)
                else:
                    AP1 = 0
                    AP1noise = 0
                    AP2 = 0
                    AP2noise = 0
                #subtract noise
                AP1subtracted = AP1 - AP1noise
                AP2subtracted = AP2 - AP2noise
    
                #write to file:
                fout.write("%d %f %f %f %f %f %f\n"%(i,AP1subtracted,AP1,
                                                   AP1noise,AP2subtracted,AP2,AP2noise))


def main():
    runs = ["20170803_1"]
    
   

    for run in runs:
        XeBCAP('/bluearc/storage/SBC-17-data/',run,'/nashome/b/bressler/sbcoutput/',1e3,120e3)
        print("run %s finished"%run)

if __name__=="__main__":
    main()