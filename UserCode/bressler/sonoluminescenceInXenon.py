#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:53:16 2020

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
import gc



def findCompression(datapath,run):
    
    runpath = datapath+"/"+run+'/'
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    
    allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
    xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
    indices = [i for i,x in enumerate(xyzf["runid"]) if str(x[0])+"_"+str(x[1]) == run]
    runposreco = {"ev":[xyzf["ev"][indices]],"x":[xyzf["bubX"][indices]],
                 "y":[xyzf["bubY"][indices]],"z":[xyzf["bubZ"][indices]]}
    pmtdiffs = []
    pmtnobubdiffs = []    
    for event in events:
        #plt.figure()
        if int(event) < 101:
            nobub = np.isnan(runposreco["z"][0][int(event)])
            e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
            sd = e["slowDAQ"]
            t = sd["elapsed_time"]
            pressure = sd["PT6"]
            trig_time = t[list(e["slowDAQ"]["TriggerOut"]).index(1.0)]
            trig_pressure = pressure[list(e["slowDAQ"]["TriggerOut"]).index(1.0)]
            
            """
            plt.figure()
            plt.plot(t,pressure)
            plt.vlines(trig_time,0,200)
            plt.show()
            """
            cgate = e["fastDAQ"]["CAMgate"]
            dcam = np.diff(cgate)
            fdt = e["fastDAQ"]["time"]
            camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]
            pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
            d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
            pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
            tracetimes = pmttracetime - pmtalign

            for t in (tracetimes):
                if t>0 and t<1:
                    lastCamOff = 0
                    for k in range(len(camOffTimes)):
                        if t+trig_time > camOffTimes[k]:
                            lastCamOff = camOffTimes[k]
                        elif t+trig_time < camOffTimes[k]:
                            break
                    if t+trig_time-lastCamOff > 25e-6:
                        if not nobub:
                            pmtdiffs.append(t)
                        elif nobub:
                            pmtnobubdiffs.append(t)
    plt.figure()
    plt.hist(pmtdiffs,50)
    plt.show()
    plt.figure()
    plt.hist(pmtnobubdiffs,50)
    plt.show()

def main():
    findCompression('/bluearc/storage/SBC-17-data/','20170708_0')

if __name__ == "__main__":
    main()