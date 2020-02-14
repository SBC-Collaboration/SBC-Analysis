#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:10:25 2020

@author: bressler
"""


import SBCcode as sbc
import numpy as np
import runlistscatalogue as rlc
import os
import pulse_integrator as pi
import matplotlib.pyplot as plt
import gc
CONVERSION_TO_CHARGE = (125.0/128)*(1/50.0)*(1/1000.0)*(1/(1.602e-19))
datapath = '/bluearc/storage/SBC-17-data/'
"""
runs = rlc.CfCombinedMultiTemp

with open("/nashome/b/bressler/sbcoutput/CfCombinedMultiTemp_merged.txt","r") as fin:
    data = fin.readlines()
    
bubdata = {} 

headers = data[0].split()

runind = headers.index("run")
evind = headers.index("event")  
xind = headers.index("x")
yind = headers.index("y")
zind = headers.index("z")
pmtt0ind = headers.index("PMTt0")
pmttraceind = headers.index("PMTt0_ind")
lagind = headers.index("lag")
spectind = headers.index("PMTphe")
blockedind = headers.index('isBlocked')
print(headers)
x=[]
y=[]
z=[]
lag = []
spect = []
goodsrcbubbles = 0
i=0
for line in data:
    if i>0:
        split_line = line.split()

        if (not np.isnan(float(split_line[xind]))) and float(split_line[blockedind])==0:
            goodsrcbubbles += 1
            evid = split_line[runind]+"-"+split_line[evind]
            #print(evid)
            x.append(float(split_line[xind]))
            y.append(float(split_line[yind]))
            z.append(float(split_line[zind]))
            lag.append(float(split_line[lagind]))
            spect.append(float(split_line[spectind]))
            bubdata[evid]=[float(split_line[spectind]),float(split_line[pmttraceind])]

    i+=1

#print(bubdata)
"""
run='20170703_2'
runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
runpath = datapath+run
events = [evnt for evnt in os.listdir(runpath) if not os.path.isfile(os.path.join(runpath,evnt))]
Nevents = len(events)
pmtfilename = runreconpath+"PMTpulseAnalysis_%s.bin"%run
pmt = sbc.DataHandling.ReadBinary.ReadBlock(pmtfilename)
print(len(pmt["PMT_pulse_area"]))
print(pmt.keys())
#evid = pmt["runid"]+"-"+pmt["ev"]
#area = pmt["PMT_pulse_area"]
#bubdata[]
evdone = []
myphe = np.zeros([len(pmt["PMT_pulse_area"]),1])
theirphe = np.zeros([len(pmt["PMT_pulse_area"]),1])

for n in range(Nevents):
    i=0
    print(n)
    indices = [i for i,x in enumerate(pmt["ev"]) if int(x)==n]
    #print(indices)
    areas=pmt["PMT_pulse_area"][indices]
    print(len(areas))
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,n)
    #print(len(e["PMTtraces"]["traces"]))
    for k in range(len(e["PMTtraces"]["traces"])):
        trace = np.fabs(e["PMTtraces"]["traces"][k][0]) 
        #if ch0 saturated, stitch in low res channel:
        if max(trace) == 128:
            trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][k][1]))
        dt = e["PMTtraces"]["dt"][k][0]
                                        
        #integrate and convert to phe:
        [a,n,totInt,pktimes] = pi.SBC_pulse_integrator_bressler(trace,dt) 
        myphe[i]=a
        theirphe[i]=-areas[i]*CONVERSION_TO_CHARGE
        i+= 1
        gc.collect()
    
plt.figure()
plt.scatter(myphe,theirphe)
plt.show()

plt.figure()
plt.hist(theirphe,30)
plt.show()
    
