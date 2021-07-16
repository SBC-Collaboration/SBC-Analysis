#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 13:16:39 2020

@author: bressler
"""
import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import pulse_integrator as pi
import runlistscatalogue as rlc
from random import randrange
from LED_blocking_scintillation import isBlocked
import gc

def randomized_t0_scintillation(runs,m):
    #m=4e7

    #m=get_gain("/bluearc/storage/SBC-17-data/",runs[0])
    
    Ncoinc = 0
    Nblocked = 0
    ntotcoinc = 0
    totevents = 0
    totbub = 0
    diffs = []
    goodz=[]
    pmtdiffs = []
    coincspec = []

    allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
    xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
    for run in runs:
        print("randomized_t0_scintillation processing run "+run)
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


        for x in events:
            totevents += 1
            if not np.isnan(runposreco["z"][0][int(x)]):
                totbub += 1
                e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,x)
                cgate = e["fastDAQ"]["CAMgate"]
                dcam = np.diff(cgate)
                fdt = e["fastDAQ"]["time"]
                camOffTimes = [fdt[i] for i in range(len(dcam)) if dcam[i] > 0.5]

                pmttracetime = e["PMTtraces"]["t0_sec"][:,0]+e["PMTtraces"]["t0_frac"][:,0]
                d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
                pmtalign = d["PMT_trigt0_sec"]+d["PMT_trigt0_frac"]
                tracetimes = pmttracetime - pmtalign
                #at0 = bubt0[int(x),0]
                rand_at0 = fdt[randrange(np.shape(fdt)[0])]
                if rand_at0 - fdt[0] < 500e-6:
                    print("too close to beginning")
                    rand_at0 = fdt[randrange(np.shape(fdt)[0])] # just trying one more time
                #print("real t0: %f"%bubt0[int(x),0])
                #print("randomized t0: %f"%rand_at0)
                i=0 # to match the indexing of the pre-made code I had 1???
                candidate = 0
                candidate_time = 0
                #candidate_PMTtime = 0
                #candidate_index = 0
                for t in (tracetimes-rand_at0):
                    # loop through every PMT trace for the event
                    
                    if t<0 and t>-500e-6: 
                        
                        lastCamOff = 0
                        for k in range(len(camOffTimes)):
                            if t+rand_at0 > camOffTimes[k]:
                                lastCamOff = camOffTimes[k]
                            elif t+rand_at0 < camOffTimes[k]:
                                break
                        if t+rand_at0-lastCamOff > 25e-6:
                            # if the trace time is within 500 microseconds before acoustic t0:
                            ntotcoinc+=1
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
                                    #candidate_PMTtime = t+at0
                                    #candidate_index = i
                    i+=1
                #i.e. if there is a candidate PMT trace with area greater than zero
                if candidate > 0:
                    Ncoinc += 1
                    #ind = candidate_index
                    #pmtt = candidate_PMTtime
                    diffs.append(candidate_time)
                    goodz.append(runposreco["z"][0][int(x)])
                    coincspec.append(candidate)
                elif candidate == 0 and not np.isnan(float(x)): 
                    bl=isBlocked(run,int(x),float(rand_at0),float(runposreco["z"][0][int(x)]))
                    if bl: 
                        Nblocked += 1
                    else:
                        coincspec.append(0)

            gc.collect()
    print("total number of events: "+str(totevents))
    print("total number of bubbles: "+str(totbub))               
    print("total coincident triggers: "+str(ntotcoinc))
    print("total blocked bubbles: "+str(Nblocked))
    print("total coincident bubbles with scintillation greater than 0phe: "+str(Ncoinc))
    print("fraction of bubbles with a coincident scintillation signal greater than 0phe: "+str(Ncoinc*100/totbub)+"%")
    
    return [coincspec,Ncoinc,ntotcoinc,totevents,totbub]


def main():
    randomized_t0_scintillation(rlc.cfJune27and28,4e7)
    
    
if __name__ == "__main__":
    main()