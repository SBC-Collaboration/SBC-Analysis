#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 09:37:12 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
from gaincalc import get_gain
import pulse_integrator as pi
from runlistscatalogue import *
import gc

CONVERSION_TO_CHARGE = (125.0/128)*(1/50.0)*(1/1000.0)*(1/(1.602e-19))

def trig_difference(runs):
    pmtdiffs = []
    pmtnobubdiffs = []
    dubbubdiffs = []
    for run in runs:
        print(run)
        runrawpath = "/bluearc/storage/SBC-17-data/%s/"%run
        runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
        acousticfilename = runreconpath+"AcousticAnalysis_%s.bin"%run
        getbubfile = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/HumanGetBub_%s.bin"%run
        a = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
        c = sbc.DataHandling.ReadBinary.ReadBlock(getbubfile)
        eventn = c["ev"]
        bubt0 = a["bubble_t0"]
        #events = [evnt for evnt in listdir(runrawpath) if not isfile(join(runrawpath,evnt))]
        #for x in events:
        for x in range(101):
            gc.collect()
            try:
                e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,x,'fastDAQ','PMTtraces')
                cgate = e["fastDAQ"]["CAMgate"]
                dcam = np.diff(cgate)
                fdt = e["fastDAQ"]["time"]
                camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]
                pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
                d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
                pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
                tracetimes = pmttracetime - pmtalign
                
                at0 = bubt0[int(x),0]
                for t in (tracetimes-at0):
                    if t<0 and t>-500e-6:
                        lastCamOff = 0
                        for k in range(len(camOffTimes)):
                            if t+at0 > camOffTimes[k]:
                                lastCamOff = camOffTimes[k]
                            elif t+at0 < camOffTimes[k]:
                                break
                        if t+at0-lastCamOff > 25e-6:
                            if list(eventn).count(int(x)) == 2:
                                pmtdiffs.append(t)
                            elif list(eventn).count(int(x)) == 1:
                                pmtnobubdiffs.append(t)
                            elif list(eventn).count(int(x)) == 3:
                                print(3)
                            elif list(eventn).count(int(x)) == 4:
                                dubbubdiffs.append(t)
            except:
                print("Last event: %d"%(x-1))
                break

    return [pmtnobubdiffs,pmtdiffs,dubbubdiffs]
    

def zdependence(runs):
    #m=4e7

    m=get_gain("/bluearc/storage/SBC-17-data/",runs[0])
    
    Ncoinc = [0,0]
    ntotcoinc = [0,0]
    totevents = [0,0]
    totbub = [0,0]
    diffs = [[],[]]
    goodz=[[],[]]
    pmtdiffs = [[],[]]
    coincspec = [[],[]]

    allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
    xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
    for run in runs:
        print("zdependence processing run "+run)
        indices = [i for i,x in enumerate(xyzf["runid"]) if str(x[0])+"_"+str(x[1]) == run]
        runposreco = {"ev":[xyzf["ev"][indices]],"x":[xyzf["bubX"][indices]],
                      "y":[xyzf["bubY"][indices]],"z":[xyzf["bubZ"][indices]]}
        runrawpath = "/bluearc/storage/SBC-17-data/%s/"%run
        runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
        acousticfilename = runreconpath+"AcousticAnalysis_%s.bin"%run
        a = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
        #c = sbc.DataHandling.ReadBinary.ReadBlock(getbubfile)
        bubt0 = a["bubble_t0"]
        events = [evnt for evnt in listdir(runrawpath) if not isfile(join(runrawpath,evnt))]
        for j in [0,1]:
            with open("/nashome/b/bressler/sbcoutput/%s_PMTmatching_ch%s.txt"%(run,str(j)),"w+") as f:
                f.write("run event PMT_t0_index PMT_t0_-at0_us PMT_t0 at0 phe z\n")
                for x in events:
                    
                    totevents[j] += 1
                    if not np.isnan(runposreco["z"][0][int(x)]):
                        totbub[j] += 1
                        e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,x)
                        veto = e["fastDAQ"]["VetoCoinc"]
                        cgate = e["fastDAQ"]["CAMgate"]
                        dcam = np.diff(cgate)
                        fdt = e["fastDAQ"]["time"]
                        camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]
    
                        
                        if min(veto)<-0.3:
                            print("Veto Coincidence")
                        pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
                        d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
                        pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
                        tracetimes = pmttracetime - pmtalign
                        at0 = bubt0[int(x),j]
                        i=1 # to match the indexing of the pre-made code 
                        candidate = 0
                        candidate_time = 0
                        candidate_PMTtime = 0
                        candidate_index = 0
                        for t in (tracetimes-at0):
                            # loop through every PMT trace for the event
                            
                            if t<0 and t>-500e-6: 
                                
                                lastCamOff = 0
                                for k in range(len(camOffTimes)):
                                    if t+at0 > camOffTimes[k]:
                                        lastCamOff = camOffTimes[k]
                                    elif t+at0 < camOffTimes[k]:
                                        break
                                if t+at0-lastCamOff > 25e-6:
                                    # if the trace time is within 500 microseconds before acoustic t0:
                                    ntotcoinc[j]+=1
                                    pmtdiffs.append(t)
                                    
                                    #take abs to get positive area:
                                    trace = np.fabs(e["PMTtraces"]["traces"][i][0]) 
                                    #if ch0 saturated, stitch in low res channel:
                                    if max(trace) == 128:
                                        trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][i][1]))
                                    dt = e["PMTtraces"]["dt"][i][0]
                                    
                                    #subtract baseline:
                                    #Actually this gets done in pulse_integrator anyway
                                    #baseline = np.mean(trace[0:50])
                                    #trace -= baseline 
                                                                
                                    #integrate and convert to phe:
                                    [phe,n,totInt,pktimes] = pi.SBC_pulse_integrator_bressler(trace,dt) 
                                    if phe != None:
                                        phe /= m
                                        #keep track of largest candidate:
                                        if phe > candidate:
                                            candidate = phe
                                            candidate_time = t
                                            candidate_PMTtime = t+at0
                                            candidate_index = i
                            i+=1
                        #i.e. if there is a candidate PMT trace with area greater than zero
                        if candidate > 0:
                            Ncoinc[j] += 1
                            ind = candidate_index
                            pmtt = candidate_PMTtime
                            diffs[j].append(candidate_time)
                            goodz[j].append(runposreco["z"][0][int(x)])
                            coincspec[j].append(candidate)
                            f.write("%s %s %d %f %f %f %f %f\n"%(run,x,ind,
                                                              candidate_time*1e6,
                                                              pmtt,at0,candidate,
                                                              runposreco["z"][0][int(x)]))
                    gc.collect()
            print("run "+run+" file %s written"%str(j))
                        #pmtdiffs.append(candidate_times[candidates.index(max(candidates))])
    print("total number of events: "+str(totevents))
    print("total number of bubbles: "+str(totbub))               
    print("total coincident triggers: "+str(ntotcoinc))
    print("total coincident bubbles with scintillation greater than 0phe: "+str(Ncoinc))
    print("fraction of bubbles with a coincident scintillation signal greater than 0phe: "+str(sum(Ncoinc)*100/sum(totbub))+"%")
    
    return [goodz,diffs,coincspec,Ncoinc,ntotcoinc,totevents,totbub]
    

    
def main():
    bgruns = bgOct10and11
    
    biberuns = BiBeSept23and24
    
    """
    ch1 files made for:
        "20171003_4","20171003_5",
        "20171004_0","20171004_1","20171004_2","20171004_3","20171004_4",
        "20171005_0","20171005_1","20171005_2","20171005_3","20171005_4",
        "20171006_0","20171006_1"
        
    ch0 files made for:
        "20171003_4","20171003_5",
        "20171004_0","20171004_1","20171004_2","20171004_3","20171004_4",
        "20171005_0","20171005_1","20171005_2","20171005_3","20171005_4",
        "20171006_0","20171006_1"
        
    """
    
    BiAlruns = []
    
    """
    ch1 files made for:
        "20171006_2","20171006_3","20171006_4","20171006_5","20171007_0",
        "20171007_1","20171007_2","20171007_4","20171007_5",
        "20171007_6","20171008_0","20171008_1","20171008_4","20171008_6",
        "20171008_7","20171009_0","20171009_1","20171009_2"
    
    ch0 files made for:
        "20171006_2","20171006_3","20171006_4","20171006_5","20171007_0",
        "20171007_1","20171007_2","20171007_4","20171007_5",
        "20171007_6","20171008_0","20171008_1","20171008_4","20171008_6",
        "20171008_7","20171009_0","20171009_1","20171009_2"
    """

    cfruns = ["20170711_15","20170711_16"]
    """
    ch0 files made for:
        "20170707_6","20170707_7","20170707_8","20170707_9","20170707_10","20170708_0",
            "20170708_1","20170708_3","20170708_5","20170708_6","20170708_7",
            "20170708_8","20170708_9","20170709_0","20170709_1","20170709_2",
            "20170709_3","20170709_4","20170709_6","20170709_7","20170709_8",
            "20170710_0","20170710_1","20170710_2","20170710_3","20170710_4",
            "20170710_5","20170710_6","20170710_7","20170710_8","20170710_9",
            "20170711_0","20170711_14","20170711_16"
    
    ch1 files made for:
        
            "20170707_6","20170707_7","20170707_8","20170707_9","20170707_10","20170708_0",
            "20170708_1","20170708_3","20170708_4","20170708_5","20170708_6",
            "20170708_7","20170708_8","20170708_9","20170709_0","20170709_1","20170709_2",
            "20170709_3","20170709_4","20170709_6",20170709_7","20170709_8",
            "20170710_0","20170710_1","20170710_2","20170710_3","20170710_4",
            "20170710_5","20170710_6","20170710_7","20170710_8", "20170710_9","20170711_0",
            "20170711_14","20170711_15","20170711_16"
            
            
    
    bad AcousticAnalysis_ .bin files:
        "20170708_2","20170708_4",
    """
    
    pmtnobubdiffs,pmtdiffs,dubbubdiffs = trig_difference(biberuns)
    goodz,diffs,coincspec,Ncoinc,ntotcoinc,totevents,totbub = zdependence(biberuns)
    
    bgpmtnobubdiffs,bgpmtdiffs,bgdubbubdiffs = trig_difference(bgruns)
    bggoodz,bgdiffs,bgcoincspec,bgNcoinc,bgntotcoinc,bgtotevents,bgtotbub = zdependence(bgruns)
    
    """
    plt.figure()
    _,bins,_=plt.hist(pmtdiffs,150,histtype='step',label="one bubble",lw=4)
    plt.hist(diffs,200,histtype='step')
    plt.hist(pmtnobubdiffs,bins,histtype='step', label = "no bubble",lw=4)
    plt.hist(dubbubdiffs,bins,histtype='step', label="two bubbles",lw=4)
    plt.xlabel("PMT trigger times minus acoustic t_0",fontsize=25)
    #plt.xlim([-100e-6,300e-6])
    plt.yscale('log')
    plt.legend(fontsize=18)
    plt.show
    """
    def ft(x,m,b):
        return m*x +b
        
    params,params_cov = scipy.optimize.curve_fit(ft,goodz,diffs)
    p1 = params[0]
    p0 = params[1]
    print("the slope is"+str(p1))
    print("The speed of sound is "+str((1/np.fabs(p1))/100)+" m/s")
    
    plt.figure()
    plt.scatter(goodz,diffs)
    plt.plot(np.arange(-3,0.1,0.1),p0+p1*np.arange(-3,0.1,0.1),lw=3)
    plt.ylim([-500e-6, 0])
    plt.xlabel("z position (cm)")
    plt.ylabel("time difference (seconds)")
    plt.show
    
    
    
    plt.figure()
    vals,bins,_=plt.hist(coincspec,np.ceil(max(coincspec)),histtype='step',color='r')
    bgvals,_,_ = plt.hist(bgcoincspec,bins,histtype='step')
    plt.xlabel("spectrum of PMT pulse areas within 500 microseconds before acoustic t_0 (photoelectrons)")
    plt.show

    plt.figure()
    plt.bar(bins[:(len(vals))],[v/totbub for v in vals],1,color='r',linewidth=0,label="californium")
    plt.bar(bins[:len(vals)],[v/bgtotbub for v in bgvals],0.7,color='b',linewidth = 0,label="background")
    plt.xlabel("spectrum of PMT pulse areas within 500 microseconds before acoustic t_0 (photoelectrons)",fontsize=25)
    plt.ylabel("probability of a scintillation pulse of this area per bubble",fontsize=25)
    plt.legend(fontsize=18)
    plt.show
    
if __name__ == "__main__":
    main()
