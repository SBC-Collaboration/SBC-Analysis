#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:58:35 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
from gaincalc import get_gain


def getRate(xyzf, runs, LTpreq):
    LT = []
    totbub = 0
    bubperpset = []
    setpoints = []
    elt = 0
    for run in runs:
        indices = [i for i,x in enumerate(xyzf["runid"]) if str(x[0])+"_"+str(x[1]) == run]
        runposreco = {"ev":[xyzf["ev"][indices]],"x":[xyzf["bubX"][indices]],
                      "y":[xyzf["bubY"][indices]],"z":[xyzf["bubZ"][indices]]}
        
        runrawpath = '/bluearc/storage/SBC-17-data/'+run
        runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
        historyfilename = runreconpath+"HistoryAnalysis_%s.bin"%run
        history = sbc.DataHandling.ReadBinary.ReadBlock(historyfilename)
    
        for eventn in range(101):
            try:
                nobub = np.isnan(runposreco["z"][0][int(eventn)])
                
                e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,eventn,"slowDAQ","event")
                pset = e["event"]["Pset"]
                elt += e["event"]["livetime"]
                if pset not in setpoints:
                    setpoints.append(pset)
                    LT.append(0)
                    bubperpset.append(0)
                ind = setpoints.index(pset)
                if not nobub:
                    totbub += 1
                    bubperpset[ind]+=1
                edges = history["PressureEdge"][0]
                times = history["PressureBins"][eventn][0][:]
                centers = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
 
                for i in range(len(centers)):
                    if abs(centers[i]-pset)<LTpreq or centers[i] < pset:
                        LT[ind]+=times[i]

            except Exception as x:
                print(x)
                break
    totlt = sum(LT)
    r = totbub/totlt
    rerr = np.sqrt(totbub)/totlt

    rperpset = np.divide(bubperpset,LT)
    r_err_perpset = np.divide(np.sqrt(bubperpset),LT)
    """
    for i in range(len(rperpset)):
        print("Pressure bin: "+str(setpoints[i]))
        #print("Live Time: "+str(LT[i]))
        #print("Bubbles: "+str(bubperpset[i]))
        print("Rate: " +str(rperpset[i])+"+/-"+str(r_err_perpset[i]))
        #print("\n")
    """
    print("25psia bin: ")
    print(str(LT[setpoints.index(25.0)])+"s")
    print(str(bubperpset[setpoints.index(25.0)])+" bubbles")
    print("rate: "+str(rperpset[setpoints.index(25.0)])+"+/-"+str(r_err_perpset[setpoints.index(25.0)])+" Hz")
    print("\n")
    print("overall rate: "+str(r)+"+/-"+str(rerr)+"Hz")
    print("total bubbles: "+str(totbub))
    print("total live time: "+str(totlt))
    print("total live time from events: " +str(elt))
    print("\n")
    return [totbub,totlt,rperpset,r_err_perpset,r,rerr,setpoints]

def main():
    allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
    xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
    
    BiBe = ["20171003_4","20171003_5","20171004_0","20171004_1","20171004_2",
            "20171004_3","20171004_4","20171005_0","20171005_1","20171005_2",
            "20171005_3","20171005_4","20171006_0","20171006_1"]
    BiAl = ["20171006_2","20171006_3","20171006_4","20171006_5",
            "20171007_0","20171007_1","20171007_2","20171007_3","20171007_4",
            "20171007_5","20171007_6","20171008_0","20171008_1","20171008_2",
            "20171008_4","20171008_5","20171008_6","20171008_7",
            "20171009_0","20171009_1","20171009_2"]
    bg = ["20171002_8","20171002_9","20171003_0"]
    YBeO = ["20170930_7","20170930_8","20171001_0","20171001_1"]
    YAl = ["20171001_6","20171001_7","20171001_8","20171002_0"]
    Co = ["20171002_5"]
    [rperpset,r_err_perpset,r,rerr,psets] = getRate(xyzf,BiAl,5)


    
if __name__ == "__main__":
    main()
