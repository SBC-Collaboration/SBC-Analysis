#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 15:31:51 2021

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
import csv


def coincidences(runs, m, PMTsearchwindow):
    
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
        bubt0 = a["bubble_t0"]
        events = [evnt for evnt in listdir(runrawpath) if not isfile(join(runrawpath,evnt))]
        for j in [0,1]:
            with open("/nashome/b/bressler/sbcoutput/%s_PMTmatching_ch%s.txt"%(run,str(j)),"w+") as f, open("/nashome/b/bressler/sbcoutput/%s_muonCoincidences.txt"%run,'w+') as fmu:
                f.write("run event PMT_t0_index PMT_t0_-at0_us PMT_t0 at0 phe z\n")
                fmu.write("run event phe\n")
                for x in events:
                    oktokeepgoing = False
                    if int(x)<len(runposreco["z"][0])-1:
                        oktokeepgoing=True
                    totevents[j] += 1
                    if oktokeepgoing and not np.isnan(runposreco["z"][0][int(x)]):
                        totbub[j] += 1
                        e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,x)
                        veto = e["fastDAQ"]["VetoCoinc"]

                        fdt = e["fastDAQ"]["time"]
                        muon = False
                        pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
                        d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
                        pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
                        tracetimes = pmttracetime - pmtalign
                        at0 = bubt0[int(x),j]
                        i=0 
                        candidate = 0
                        candidate_time = 0
                        candidate_PMTtime = 0
                        candidate_index = 0
                        for t in (tracetimes-at0):
                            # loop through every PMT trace for the event
                            
                            if t<PMTsearchwindow[1] and t>PMTsearchwindow[0]: 
                                # if the trace time is in the window
                                if max(veto)>0.1:
                                    if fdt[list(veto).index(max(veto))]-at0<PMTsearchwindow[1] and fdt[list(veto).index(max(veto))]-at0>PMTsearchwindow[0]:
                                        print("Veto Coincidence: event "+run+"-"+str(x))
                                        muon = True
                                        

                                ntotcoinc[j]+=1
                                pmtdiffs.append(t)
                                
                                #take abs to get positive area:
                                trace = np.fabs(e["PMTtraces"]["traces"][i][0]) 
                                #if ch0 saturated, stitch in low res channel:
                                if max(trace) == 128:
                                    trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][i][1]))
                                dt = e["PMTtraces"]["dt"][i][0]
                                                            
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
                        if muon:
                            fmu.write("%s %s %f\n"%(run,x,candidate))
                    gc.collect()
            print("run "+run+" file %s written"%str(j))
    print("total number of events: "+str(totevents))
    print("total number of bubbles: "+str(totbub))               
    print("total coincident triggers: "+str(ntotcoinc))
    print("total coincident bubbles with scintillation greater than 0phe: "+str(Ncoinc))
    print("fraction of bubbles with a coincident scintillation signal greater than 0phe: "+str(sum(Ncoinc)*100/sum(totbub))+"%")
    
def coincidences_handscanned(runs, m, datasetname, PMTsearchwindow):
    
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
        
        bubt0_dict = {}
       
        with open("/nashome/b/bressler/sbcoutput/handscans_acoustic_t0/acoustict0handscan_%s.csv"%datasetname, 'r') as csvfile:
            r = csv.reader(csvfile)
            i=0
            for row in r:
                if i>0:
                    bubt0_dict[row[0]+'-'+row[1]] = float(row[2])
                i+=1

        events = [evnt for evnt in listdir(runrawpath) if not isfile(join(runrawpath,evnt))]
        j=1
        with open("/nashome/b/bressler/sbcoutput/%s_PMTmatching_ch%s.txt"%(run,str(j)),"w+") as f, open("/nashome/b/bressler/sbcoutput/%s_muonCoincidences.txt"%run,'w+') as fmu:
            f.write("run event PMT_t0_index PMT_t0_-at0_us PMT_t0 at0 phe z\n")
            fmu.write("run event phe\n")
            
            for x in events:
                okaytokeepgoing = False
                if int(x)<len(runposreco["z"][0])-1:
                    okaytokeepgoing=True
                totevents[j] += 1
                if okaytokeepgoing and not np.isnan(runposreco["z"][0][int(x)]):
                    totbub[j] += 1
                    e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,x)
                    veto = e["fastDAQ"]["VetoCoinc"]

                    fdt = e["fastDAQ"]["time"]
                    muon = False
                    pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
                    d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
                    pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
                    tracetimes = pmttracetime - pmtalign
                    at0 = bubt0_dict[run+'-'+str(x)]
                    i=0 
                    candidate = 0
                    candidate_time = 0
                    candidate_PMTtime = 0
                    candidate_index = 0
                    for t in (tracetimes-at0):
                        # loop through every PMT trace for the event
                        
                        if t<PMTsearchwindow[1] and t>PMTsearchwindow[0]: 
                            # if the trace time is in the window
                            if max(veto)>0.1:
                                if fdt[list(veto).index(max(veto))]-at0<PMTsearchwindow[1] and fdt[list(veto).index(max(veto))]-at0>PMTsearchwindow[0]:
                                    print("Veto Coincidence: event "+run+"-"+str(x))
                                    muon = True
                                    

                            ntotcoinc[j]+=1
                            pmtdiffs.append(t)
                            
                            #take abs to get positive area:
                            trace = np.fabs(e["PMTtraces"]["traces"][i][0]) 
                            #if ch0 saturated, stitch in low res channel:
                            if max(trace) == 128:
                                trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][i][1]))
                            dt = e["PMTtraces"]["dt"][i][0]
                                                        
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
                    if muon:
                        fmu.write("%s %s %f\n"%(run,x,candidate))
                gc.collect()
            print("run "+run+" file %s written"%str(j))
    print("total number of events: "+str(totevents))
    print("total number of bubbles: "+str(totbub))               
    print("total coincident triggers: "+str(ntotcoinc))
    print("total coincident bubbles with scintillation greater than 0phe: "+str(Ncoinc))
    print("fraction of bubbles with a coincident scintillation signal greater than 0phe: "+str(sum(Ncoinc)*100/sum(totbub))+"%")
    
    
def find_coincidence(runs, m, isHandscannedT0, datasetname, PMTsearchwindow):
    if isHandscannedT0:
        coincidences_handscanned(runs, m, datasetname, PMTsearchwindow)
    else:
        coincidences(runs, m, PMTsearchwindow)