#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 13:47:20 2020

@author: bressler
"""
import SBCcode as sbc
import matplotlib.pyplot as plt
from PICOcode.REFPROP.SeitzModel import SeitzModel

def plot_xebc_slowDAQ(runpath,event):
    indices_back = 30
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event,"slowDAQ","event")
    sd = e["slowDAQ"]
    t =sd["elapsed_time"]
    T1 = sd["T1"]
    PT4 = sd["PT4"]
    PT5 = sd["PT5"]
    PT6 = sd["PT6"]
    ev = e["event"]
    pset = ev["Pset"]
    trig_slow = ev["trigger_slowDAQ"]
    trig_cam =ev["trigger_cameras"]
    trig_plc = ev["trigger_PLC"]
    trig_pressure = e["slowDAQ"]["PT6"][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
    fig=plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    plt.title(runpath.split('/')[4])
    ax1.plot(t,PT4,label="PT4")
    ax1.plot(t,PT5,label="PT5")
    ax1.plot(t,PT6,label="PT6")
    ax1.plot([t[0],t[-1]],[pset,pset],label="Pressure Setpoint")
    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('Pressure [psia]')
    ax1.legend(fontsize=18)
    ax1.text(5,50,"slowDAQ: "+str(trig_slow))
    ax1.text(5,80,"cam: "+str(trig_cam))
    ax1.text(5,110,"PLC: "+str(trig_plc))
    
    ax2 = fig.add_subplot(1,2,2)
    ax2.plot(t,T1,label='T1')
    ax2.legend(fontsize=18)
    ax2.set_xlabel('time [s]')
    ax2.set_ylabel('Temperature [C]')
    plt.show()
    
def main():
    plot_xebc_slowDAQ('/bluearc/storage/SBC-17-data/20170629_7',0)

if __name__=="__main__":
    main()
    
    