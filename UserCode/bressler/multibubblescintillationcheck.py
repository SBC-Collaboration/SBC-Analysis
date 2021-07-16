#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 19:33:02 2021

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import pulse_integrator as pi
import gc

def check_multibub_scintillation(run, event, at0, PMTgain, PMTwindow):
    tstart = PMTwindow[0]
    t_end= PMTwindow[1]
    scintillation_signal = 0
    datadir = '/bluearc/storage/SBC-17-data'
    e = sbc.DataHandling.GetSBCEvent.GetEvent(datadir+'/'+run,event)
    cgate = e["fastDAQ"]["CAMgate"]
    fdt = e["fastDAQ"]["time"]
    
    LED_on = [fdt[i] for i in range(len(cgate)) if cgate[i]<-0.5]
    look_times = [x for x in LED_on if (x < 0 and abs(x-at0)<tstart)]
    #print(str(len(LED_on)/len(fdt)))
    if len(look_times)>0:
        LED_during_window = True
        
    else:
        LED_during_window = False

        dcam = np.diff(cgate)
        fdt = e["fastDAQ"]["time"]
        camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]
        pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
        d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
        pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
        tracetimes = pmttracetime - pmtalign
        
        i=0 # to match the indexing of the pre-made code I had 1???
        candidate = 0

        for t in (tracetimes-at0):
            # loop through every PMT trace for the event
            
            if t<t_end and t>tstart: 
                # if the trace time is within 500 microsec before acoustic t0
                """        
                lastCamOff = 0
                for k in range(len(camOffTimes)):
                    if t+at0 > camOffTimes[k]:
                        lastCamOff = camOffTimes[k]
                    elif t+at0 < camOffTimes[k]:
                        break
                if t+at0-lastCamOff > 25e-6:
                    # if the trace time is more than 25 microseconds away from a camera gate rise
                """    
                #take abs to get positive area:
                trace = np.fabs(e["PMTtraces"]["traces"][i][0]) 
                #if ch0 saturated, stitch in low res channel:
                if max(trace) == 128:
                    trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][i][1]))
                dt = e["PMTtraces"]["dt"][i][0]

                                            
                #integrate and convert to phe:
                [phe,n,totInt,pktimes] = pi.SBC_pulse_integrator_bressler(trace,dt) 
                if phe != None:
                    phe /= PMTgain
                    #keep track of largest candidate:
                    if phe > candidate:
                        candidate = phe
                            
            i+=1
        #i.e. if there is a candidate PMT trace with area greater than zero
        if candidate > 0:
            scintillation_signal = candidate

    gc.collect()
    return [LED_during_window, scintillation_signal]

def main():
    returned = check_multibub_scintillation('/bluearc/storage/SBC-17-data/20170703_5', 5, -0.1)
    
if __name__ == "__main__":
    main()