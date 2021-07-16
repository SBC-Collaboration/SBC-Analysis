#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:58:35 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import scipy
from gaincalc import get_gain
from PICOcode.REFPROP.SeitzModel import SeitzModel
import runlistscatalogue as rlc

OKLONGEVENTS = ['20170701_6-81', '20170627_3-55', '20170628_5-29', '20170623_7-38',
                '20170623_7-76', '20170624_0-35', '20170624_0-72', '20170801_2-39',
                '20170717_0-67', '20170626_5-6', '20170626_5-20', '20170626_9-2',
                '20170626_9-12', '20170626_9-77', '20170626_9-90', '20170627_0-31',
                '20170627_0-73', '20170627_0-74', '20170702_4-19', '20170702_4-28',
                '20170703_0-59', '20170703_6-37', '20170704_1-97']


def getRate(xyzf, runs, LTpreq,PT):
    LT = []
    totbub = 0
    bubperpset = []
    doublesperpset = []
    triplesperpset = []
    setpoints = []
    elt = 0
    feedbackTransducer = int(PT[-1])-1
    temperatures = []
    counts = []
    
    event_livetimes = []
    event_counter = 0
    for run in runs:
        runpath = '/bluearc/storage/SBC-17-data/'+run+'/'
        events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
        Nevents = len(events)
        print("run "+run+": %d events"%Nevents)
        indices = [i for i,x in enumerate(xyzf["runid"]) if str(x[0])+"_"+str(x[1]) == run]
        runposreco = {"ev":[xyzf["ev"][indices]],"x":[xyzf["bubX"][indices]],
                      "y":[xyzf["bubY"][indices]],"z":[xyzf["bubZ"][indices]]}
        
        runrawpath = '/bluearc/storage/SBC-17-data/'+run
        runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
        historyfilename = runreconpath+"HistoryAnalysis_%s.bin"%run
        history = sbc.DataHandling.ReadBinary.ReadBlock(historyfilename)
        
        getbubfile = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/HumanGetBub_%s.bin"%run
        c = sbc.DataHandling.ReadBinary.ReadBlock(getbubfile)
        #print(c.keys())
        evn = c["ev"]
        #print(evn)
        #print(c["nbubimage"])
        count = Counter(evn)
        #print(count)
        nbubs = np.zeros(Nevents)
        #print(len(nbubs))
        for c in count:
            #print(c)
            N = count[c]
            #print(N)
            counts.append(N)
            if N%2 != 0 and N != 1:
                print("lines isn't even")
                print(c)
                print(N)
            elif N%2 == 0:
                n = N/2
                nbubs[c]=n
        e0 = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,0,"slowDAQ","event")
        T25 = np.mean(e0["slowDAQ"]["T1"]) 
        seitz25 = SeitzModel(25,T25,'xenon')
        Q25 = seitz25.Q
        #print(nbubs)
        for eventn in range(Nevents):
            try:
                nobub = np.isnan(runposreco["z"][0][int(eventn)])
                n = nbubs[int(eventn)]
                #print("n=%d"%n)
                
                if nobub and n==1:
                    print("nobub and n=1: " + run +'-'+str(eventn))
                """
                elif (not nobub) and n==1:
                    print("(not nobub) and n = 1: " + run + '-'+str(eventn))
                elif n != 1 and nobub:
                    print("nobub and n != 1: "+run+'-'+str(eventn))
                elif n !=1 and not nobub:
                    print("not nobub and n != 1: "+run+'-'+str(eventn))
                """
                #print(n==1)
                e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,eventn,"slowDAQ","event")
                lt_event = e["event"]["livetime"]
                if(lt_event>800):
                    if run+'-'+str(eventn) not in OKLONGEVENTS:
                        print("check long event, excluding now: "+run+"-"+str(eventn))
                        continue
                event_livetimes.append(lt_event)
                event_counter += 1
                if lt_event < 20:
                    #print("short event excluded")
                    continue
                T = np.mean(e["slowDAQ"]["T1"])
                #print(T)
                if T>-1000:
                    temperatures.append(T)
                else:
                    print("Not including unpysical temperature %f C"%T)
                pset = e["event"]["Pset"]
                #seitz = SeitzModel(pset,T,'xenon')
                #Q = seitz.Q
                elt += e["event"]["livetime"]
                if pset not in setpoints:
                    setpoints.append(pset)
                    LT.append(0)
                    bubperpset.append(0)
                    doublesperpset.append(0)
                    triplesperpset.append(0)
                ind = setpoints.index(pset)
                
                indices_back = 30
                trig_pressure = e["slowDAQ"][PT][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
                
                #if not nobub:
                    #print("adding bubble at P = %f"%trig_pressure)
                if n > 0:
                    totbub += 1
                if (np.abs(trig_pressure-pset)<LTpreq/2):
                    if n == 1: bubperpset[ind]+=1
                    elif n==2: doublesperpset[ind] += 1
                    elif n==3: triplesperpset[ind] += 1
                    
                edges = history["PressureEdge"][feedbackTransducer]
                times = history["PressureBins"][eventn][feedbackTransducer][:]
                centers = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
                #print(edges)
                for i in range(len(centers)):
                    if abs(centers[i]-pset)<LTpreq/2:
                        LT[ind]+=times[i]

            except Exception as x:
                print(x)
                break
    totlt = sum(LT)
    r = totbub/totlt
    rerr = np.sqrt(totbub)/totlt

    rperpset = np.divide(bubperpset,LT)
    r_err_perpset = np.divide(np.sqrt(bubperpset),LT)

    for i in range(len(rperpset)):
        if LT[i] > 50:
            seitz = SeitzModel(setpoints[i],np.mean(temperatures),'xenon')
            #print(setpoints[i])
            #print(min(temperatures))
            #print(temperatures)
            #print(seitz)
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print("Pressure bin: "+str(setpoints[i]))
            print("Live Time: "+str(LT[i]))
            print("Single Bubbles: "+str(bubperpset[i]))
            print("Singles Rate: " +str(rperpset[i])+"+/-"+str(r_err_perpset[i]))
            print("Double Bubbles: %d"%doublesperpset[i])
            print("Doubles Rate: %.2e"%(doublesperpset[i]/LT[i]))
            print("Triple Bubbles: %d"%triplesperpset[i])
            print("Triples Rate: %.2e"%(triplesperpset[i]/LT[i]))
            print("Seitz Threshold: %f keV"%seitz.Q)
            print("Temperature: %.2fC"%np.mean(temperatures))
            print("\n")
    if 25.0 in setpoints:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("25 psia bin: ")
        print(str(LT[setpoints.index(25.0)])+"s")
        print(str(bubperpset[setpoints.index(25.0)])+" bubbles")
        print(str(Q25) + " keV")
        print(str(T25) + " Celsius")
        print("rate: "+str(rperpset[setpoints.index(25.0)])+"+/-"+str(r_err_perpset[setpoints.index(25.0)])+" Hz")
        print("\n")
        print("overall rate: "+str(r)+"+/-"+str(rerr)+"Hz")
        print("total events: %d"%event_counter)
        print("total bubbles: "+str(totbub))
        print("fraction of events with bubbles: %f"%(totbub/event_counter))
        print("total live time: "+str(totlt))
        print("total live time from events: " +str(elt))
        print("mean live time: %f sec"%np.mean(event_livetimes))
        print("trigger rate: %f Hz"%(1/float(np.mean(event_livetimes))))
        print("bubble rate from bub fraction  and LT: %f Hz"%(float((totbub/event_counter))/float(np.mean(event_livetimes))))
        print("\n")
    plt.figure()
    plt.scatter(np.arange(event_counter),event_livetimes, facecolor='none', edgecolor='b')
    plt.grid()
    plt.xlabel('event count')
    plt.ylabel('event livetime (s)')
    plt.show()
    plt.figure()
    plt.hist(event_livetimes,int(np.ceil(np.sqrt(len(event_livetimes)))))
    plt.show()
    return [totbub,totlt,rperpset,r_err_perpset,r,rerr,setpoints]

def main():
    allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
    xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
    
    [totbub,totlt,rperpset,r_err_perpset,r,rerr,setpoints] = getRate(xyzf,rlc.bgJuly3and4,1,"PT6")
    print(totbub)
    print(totlt)
    print(r)
    print(rerr)

    
if __name__ == "__main__":
    main()
