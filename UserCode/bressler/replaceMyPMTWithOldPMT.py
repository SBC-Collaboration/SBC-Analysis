#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:28:58 2020

@author: bressler
"""
import SBCcode as sbc
import runlistscatalogue as rlc


runsName = 'bgJuly11and12'
runs = rlc.bgJuly11and12

with open('/nashome/b/bressler/sbcoutput/%s_merged_oldPMT.txt'%runsName,'w') as fout:
    with open("/nashome/b/bressler/sbcoutput/%s_merged.txt"%runsName,"r") as fin:
            data = fin.readlines()
    
   
    headers = data[0].split()
    runind = headers.index("run")
    eventind = headers.index("event")
    xind = headers.index("x")
    yind = headers.index("y")
    zind = headers.index("z")
    lagind = headers.index("lag")
    spectind = headers.index("PMTphe")
    blockedind = headers.index('isBlocked')
    x=[]
    y=[]
    z=[]
    evid = []
    lag = []
    spect = []
    isblocked = []
    i=0
    for line in data:
        if i>0:
            split_line = line.split()
            run = split_line[runind]
            evn = int(split_line[eventind])
            runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
            timingname = runreconpath+"TimingAnalysis_%s.bin"%run
            timing = sbc.DataHandling.ReadBinary.ReadBlock(timingname)
         
            #spect.append(float(split_line[spectind]))
            print(str(timing["PMTmatch_nphe"][evn][0]))
            print("mine: %s"%str(split_line[spectind]))
            split_line[spectind] = str(timing["PMTmatch_nphe"][evn][0])
            #theirspect.append(timing["PMTmatch_nphe"][i-1][0])
            line = ("%s "*len(split_line)+"\n")%tuple(split_line)
            print(line)
            fout.write(line)
        else: 
            fout.write(line)
        i+=1