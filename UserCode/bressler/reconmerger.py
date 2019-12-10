#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 09:19:54 2019

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import math
from LED_blocking_scintillation import isBlocked
from runlistscatalogue import *

def merge(runsToMerge,name):

    allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
    xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
    totevents = 0
    blockedEvents = 0
    nosigNotBlockedEvents = 0
    with open("/nashome/b/bressler/sbcoutput/%s_merged.txt"%name,"w+") as fout,\
        open("/nashome/b/bressler/sbcoutput/%s_nosig_merged.txt"%name,"w+") as fnosig:
        fout.write("run event x y z AP1 AP2 at0_0 at0_1 PMTt0 PMTt0_ind PMTphe lag which_acoustic_channel isBlocked\n")
        
        for run in runsToMerge:
            ch0file=open("/nashome/b/bressler/sbcoutput/%s_PMTmatching_ch0.txt"%run,"r")
            ch1file=open("/nashome/b/bressler/sbcoutput/%s_PMTmatching_ch1.txt"%run,"r")
            
            APfile = open("/nashome/b/bressler/sbcoutput/XeBCAP_%s.txt"%run,"r")
            
            indices = [i for i,x in enumerate(xyzf["runid"]) if str(x[0])+"_"+str(x[1]) == run]
            runposreco = {"ev":[xyzf["ev"][indices]],"x":[xyzf["bubX"][indices]],
                          "y":[xyzf["bubY"][indices]],"z":[xyzf["bubZ"][indices]]}
            runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
            acousticfilename = runreconpath+"AcousticAnalysis_%s.bin"%run
            a = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
            bubt0 = a["bubble_t0"]
            
            AP1s = []
            AP2s = []
            
            ev_with_scintillation_0 = []
            at0_0s = []
            indices_0 = []
            pmtt0s_0 = []
            lags_0 = []
            phes_0 = []
            
            ev_with_scintillation_1 = []
            at0_1s = []
            indices_1 = []
            pmtt0s_1 = []
            lags_1 = []
            phes_1 = []
            
            ch0 = enumerate(ch0file)
            ch1 = enumerate(ch1file)
            for i,l in ch0:
                if i>=1:
                    d=l.split()
                    run_num = d[0]
                    event_num = d[1]
                    ev_with_scintillation_0.append(int(event_num))
                    ind = d[2]
                    indices_0.append(ind)
                    lag = d[3]
                    lags_0.append(lag)
                    pmtt0 = d[4]
                    pmtt0s_0.append(pmtt0)
                    at0 = d[5]
                    at0_0s.append(at0)
                    phe = d[6]
                    phes_0.append(phe)
                    z = d[7]
                    
            
            for i,l in ch1:
                if i>=1:
                    d=l.split()
                    run_num = d[0]
                    event_num = d[1]
                    ev_with_scintillation_1.append(int(event_num))
                    ind = d[2]
                    indices_1.append(ind)
                    lag = d[3]
                    lags_1.append(lag)
                    pmtt0 = d[4]
                    pmtt0s_1.append(pmtt0)
                    at0 = d[5]
                    at0_1s.append(at0)
                    phe = d[6]
                    phes_1.append(phe)
                    z = d[7]
    
        
            ch0file.close()
            ch1file.close()
            
            ap = enumerate(APfile)
            for i,l in ap:
                if i>=1:
                    d=l.split()
                    AP1s.append(float(d[1])) # 2 for raw
                    AP2s.append(float(d[4])) # 5 for raw
            
            for evnt in runposreco["ev"][0]:
                totevents += 1
                AP1 = AP1s[evnt]
                AP2 = AP2s[evnt]
                x=runposreco["x"][0][evnt]
                y=runposreco["y"][0][evnt]
                z=runposreco["z"][0][evnt]
                at0_0=bubt0[evnt,0]
                at0_1=bubt0[evnt,1]
                in0=False
                in1=False
                PMTt0 = np.nan
                sig = np.nan
                pmtt0ind = np.nan
                lag = np.nan
                chan = np.nan
                
                if evnt in ev_with_scintillation_0:
                    in0=True
                    ind0 = ev_with_scintillation_0.index(evnt)
    
                if evnt in ev_with_scintillation_1:
                    in1=True
                    ind1 = ev_with_scintillation_1.index(evnt)
                    
                    
                if (in0 and in1) or (in0 and not in1):
                    PMTt0 = pmtt0s_0[ind0]
                    sig = phes_0[ind0]
                    pmtt0ind = indices_0[ind0]
                    lag = lags_0[ind0]
                    chan = 0.
                    
                elif (in1 and not in0):
                    PMTt0 = pmtt0s_1[ind1]
                    sig = phes_1[ind1]
                    pmtt0ind = indices_1[ind1]
                    lag = lags_1[ind1]
                    chan = 1.
                bl = 0
                if math.isnan(float(sig)) and not math.isnan(float(x)):    
                    bl=isBlocked(run,int(evnt),float(at0_0),float(z))
                    if not bl:
                        fnosig.write("%s %d\n"%(str(run),int(evnt)))
                        nosigNotBlockedEvents += 1
                    else:
                        blockedEvents +=1
                        
                fout.write("%s %d %f %f %f %f %f %f %f %f %.0f %f %f %.0f %0f\n"%(str(run),int(evnt),
                    float(x),float(y),float(z),AP1,AP2,float(at0_0),float(at0_1),float(PMTt0),
                                                                  float(pmtt0ind),float(sig),
                                                                  float(lag),float(chan),float(bl)))
    print("fraction of events blocked by LED: %f"%(blockedEvents/totevents))
    print("fraction of events with no scintillation and not blocked: %f"%(nosigNotBlockedEvents/totevents))
            
def main():
    merge(BiBeOct3to6,'BiBeOct3to6')

if __name__=="__main__":
    main()