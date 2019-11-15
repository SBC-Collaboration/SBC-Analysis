#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 07:47:36 2019

@author: bressler
"""


import SBCcode as sbc
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from runlistscatalogue import *


def mult(runs):
    counts = []
    nbubs = []
    
    eventcount = np.zeros(87)
    LT = np.zeros(87)
    expand_times = [[] for i in range(87)]
    totbub = 0
    onebubcount = np.zeros(87)
    twobubcount = np.zeros(87)
    threebubcount = np.zeros(87)
    setpoints = []
    elt = 0
    
    allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
    
    for run in runs:
        print(run)
        tcut = 0
        
        runrawpath = '/bluearc/storage/SBC-17-data/'+run
        runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
        historyfilename = runreconpath+"HistoryAnalysis_%s.bin"%run
        history = sbc.DataHandling.ReadBinary.ReadBlock(historyfilename)
        getbubfile = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/HumanGetBub_%s.bin"%run
        c = sbc.DataHandling.ReadBinary.ReadBlock(getbubfile)
        
        eventn = c["ev"]
        count = Counter(eventn)
        
        edges = history["PressureEdge"][5]
        centers = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]

        for c in count:
            N = count[c]
            counts.append(N)
            if N%2 != 0 and N != 1:
                print("lines isn't even")
                print(c)
                print(N)
            elif N%2 == 0:
                n = N/2
                nbubs.append(n)
        
        for eventn in range(101):
            try:
                #onebub = not np.isnan(runposreco["z"][0][int(eventn)])
                n = 0
                N = count[eventn]
                if N%2 == 0:
                    n = N/2
                e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,eventn,"slowDAQ","event")
                t = e["slowDAQ"]["elapsed_time"]
                tcenters = [(t[i+1] + t[i])/2 for i in range(len(t)-1)]
                trig_time = t[list(e["slowDAQ"]["TriggerOut"]).index(1.0)-20]
                pslope = np.diff(e["slowDAQ"]["PT6"])
                expstartind = list(pslope).index(min(list(pslope)))
                exp_start = tcenters[expstartind]
                trig_pressure = e["slowDAQ"]["PT6"][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-20]
                pset = e["event"]["Pset"]
                ev_lt = e["event"]['livetime']
                elt += ev_lt
                if pset not in setpoints:
                    setpoints.append(pset)
                #print(len(history["PressureBins"][eventn]))
                times = history["PressureBins"][eventn][5][:]
                """
                plt.figure()
                plt.plot([trig_time,trig_time],[0,200])
                plt.plot([exp_start,exp_start],[0,200])
                plt.plot([0,max(t)],[trig_pressure,trig_pressure])
                plt.plot(t,e["slowDAQ"]["PT6"])
                plt.plot(tcenters,pslope)
                plt.xlabel("time",fontsize=18)
                plt.ylabel("signal",fontsize=18)
                plt.show
                """
                if n == 1:
                    totbub += 1
                    for i in range(len(centers)):
                        if trig_time > exp_start + tcut:
                            if trig_pressure >= edges[i] and trig_pressure <= edges[i+1]:
                                onebubcount[i] += 1
                                eventcount[i] += 1
                                
                                expand_times[i].append(ev_lt-tcut)
                            LT[i] += times[i]
                        
                        
                elif n == 2:
                    totbub += 1
                    for i in range(len(centers)):
                        if trig_time > exp_start + tcut:
                            if trig_pressure >= edges[i] and trig_pressure <= edges[i+1]:
                                eventcount[i] += 1
                                twobubcount[i] += 1
                                expand_times[i].append(ev_lt)
                            LT[i] += times[i]

                elif n >= 3:
                    for i in range(len(centers)):
                        if trig_time > exp_start + tcut:
                            if trig_pressure >= edges[i] and trig_pressure <= edges[i+1]:
                                eventcount[i] += 1
                                threebubcount[i] += 1
                                expand_times[i].append(ev_lt)
                            LT[i] += times[i]

                else:
                    for i in range(len(centers)):
                        if trig_time > exp_start + tcut:
                            if trig_pressure >= edges[i] and trig_pressure <= edges[i+1]:
                                eventcount[i] += 1
                                expand_times[i].append(ev_lt)
                            LT[i] += times[i]

            except Exception as x:
                print(x)
                break

    rateList1 = [onebubcount[i]/LT[i] for i in range(len(onebubcount))]
    rateErrList1h = [(0.5+np.sqrt(onebubcount[i]+0.25))/LT[i] for i in range(len(onebubcount))]
    rateErrList1l = [(-0.5+np.sqrt(onebubcount[i]+0.25))/LT[i] for i in range(len(onebubcount))]


    rateList2 = [twobubcount[i]/LT[i] for i in range(len(twobubcount))]
    rateErrList2h = [(0.5+np.sqrt(twobubcount[i]+0.25))/LT[i] for i in range(len(twobubcount))]
    rateErrList2l = [(-0.5+np.sqrt(twobubcount[i]+0.25))/LT[i] for i in range(len(twobubcount))]
    
    rateList3 = [threebubcount[i]/LT[i] for i in range(len(threebubcount))]
    rateErrList3h = [(0.5+np.sqrt(threebubcount[i]+0.25))/LT[i] for i in range(len(threebubcount))]
    rateErrList3l = [(-0.5+np.sqrt(threebubcount[i]+0.25))/LT[i] for i in range(len(threebubcount))]
    
    rateListE = [eventcount[i]/LT[i] for i in range(len(eventcount))]
    rateErrListEh = [(0.5+np.sqrt(eventcount[i]+0.25))/LT[i] for i in range(len(eventcount))]
    rateErrListEl = [(-0.5+np.sqrt(eventcount[i]+0.25))/LT[i] for i in range(len(eventcount))]
    
    plt.figure()
    plt.errorbar(centers,rateList1,[rateErrList1l,rateErrList1h],fmt='ro',label='Single Bubbles')
    plt.errorbar(centers,rateList2,[rateErrList2l,rateErrList2h],fmt='go',label='Double Bubbles')
    plt.errorbar(centers,rateList3,[rateErrList3l,rateErrList3h],fmt='bo',label='Triple Bubbles')
    plt.errorbar(centers,rateListE,[rateErrListEl,rateErrListEh],fmt='ko',label='All Events')
    plt.yscale('log')
    plt.xlabel('Pressure [psia]',fontsize=18)
    plt.ylabel('Rate [Hz]', fontsize=18)
    plt.legend(fontsize=18)
    plt.grid()
    plt.show
    
    plt.figure()
    hy = []
    hx = []
    for i in range(len(expand_times)):
        if len(expand_times[i])>0:
            for j in range(len(expand_times[i])):
                hx.append(centers[i])
                hy.append(expand_times[i][j])

    plt.hist2d(hx,hy,bins=(len(centers),30))
    plt.xlabel('PT6 [psia]',fontsize=18)
    plt.ylabel('Total time since expansion [s]',fontsize=18)
    plt.show
    
    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel('pressure (psia)',fontsize=18)
    ax1.set_ylabel('Live time (s)',fontsize=18)
    ax2.set_ylabel('Number of bubbles',fontsize=18)
    ax1.scatter(centers,LT,20,'r')
    ax2.scatter(centers,onebubcount,20,'b')
    plt.show
    
    totbub = len(nbubs)
    print("total number of bubbles: %d"%totbub)
    print("total number of events: %d" %len(counts))
    m = max(nbubs)
    for i in range(len(nbubs)):
        if nbubs[i]>3:
            nbubs[i] = 3
            
    binedges = [0.5,1.5,2.5,3.5]
    bincenters = [1,2,3]
    plt.figure()
    ns, _ = np.histogram(nbubs,binedges)
    plt.errorbar(bincenters,ns,np.sqrt(ns),fmt='o')
    plt.yscale('log')
    plt.title("Max Nbub: %d"%int(m),fontsize=20)
    plt.xlabel("Nbub",fontsize=18)
    plt.ylabel("counts",fontsize=18)
    plt.xlim(0.5,3.5)
    plt.grid()
    plt.show
    
    f = [ns[i]/totbub for i in range(len(ns))]
    ferrh = [(0.5+np.sqrt(ns[i]+0.25))/totbub for i in range(len(ns))]
    ferrl = [(-0.5+np.sqrt(ns[i]+0.25))/totbub for i in range(len(ns))]

    plt.figure()
    plt.errorbar(bincenters,f,[ferrl,ferrh],fmt='o')
    plt.xlabel("Nbub",fontsize=18)
    plt.ylabel("fraction of bubbles",fontsize=18)
    plt.xlim([0.5,3.5])
    plt.grid()
    plt.show
        

def main():
    runs = cfJune28to30
    mult(runs)
if __name__=="__main__":
    main()