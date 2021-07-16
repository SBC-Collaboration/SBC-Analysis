#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:35:21 2019

@author: bressler
"""


import SBCcode as sbc
import numpy as np


def isBlocked(run, event, at0, z, PMTwindow):
    #print("testing")
    ret=False
    #lag_expected = (-23.387649*z - 261.020495)*1e-6 # fit from other analysis
    #t0_expected = at0 + lag_expected
    datadir = '/bluearc/storage/SBC-17-data'
    e = sbc.DataHandling.GetSBCEvent.GetEvent(datadir+'/'+run,event)
    #p1 = e["fastDAQ"]["Piezo1"]
    cgate = e["fastDAQ"]["CAMgate"]
    #dcam = np.diff(cgate)
    fdt = e["fastDAQ"]["time"]
    #tstep = fdt[1]-fdt[0]
    #camOnTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] < -0.3]
    #camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.3]
    
    LED_on = [fdt[i] for i in range(len(cgate)) if cgate[i]<-0.5]
    look_times = [x for x in LED_on if (x < at0 + PMTwindow[1] and x > at0 + PMTwindow[0])]
    #if len(look_times) == 0:
        #print("look_times empty")
    """
    #print(str(len(LED_on)/len(fdt)))

    for t in look_times:
        if np.absolute(t-t0_expected)<25e-6:
            #print("LED blocked it")
            ret = True
            break
    """
    if len(look_times) > 0:
        ret = True
    return ret
    
def main():
    
    print(isBlocked('20171007_5',2,-0.025,-1.5))
    
if __name__=="__main__":
    main()

    