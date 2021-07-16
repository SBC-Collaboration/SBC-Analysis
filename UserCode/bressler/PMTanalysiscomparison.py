#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:10:25 2020

@author: bressler
"""


import SBCcode as sbc
import numpy as np
import scipy
import runlistscatalogue as rlc
import os
import pulse_integrator as pi
import matplotlib.pyplot as plt
plt.style.use('seaborn-darkgrid')
import gc
from gaincalc import get_gain

CONVERSION_TO_CHARGE = (125.0/128)*(1/50.0)*(1/1000.0)*(1/(1.602e-19))
datapath = '/bluearc/storage/SBC-17-data/'
"""
run='20170628_9'
runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
timingname = runreconpath+"TimingAnalysis_%s.bin"%run
timing = sbc.DataHandling.ReadBinary.ReadBlock(timingname)
print(timing.keys())

#runs = rlc.CfCombinedMultiTemp

with open("/nashome/b/bressler/sbcoutput/20170628_9_merged.txt","r") as fin:
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
theirspect = []
goodsrcbubbles = 0
i=0
for line in data:
    if i>0:
        split_line = line.split()
        if (not np.isnan(float(split_line[xind]))) and float(split_line[blockedind])==0:
            goodsrcbubbles += 1
            evid = split_line[runind]+"-"+split_line[evind]
            print(evid)
            print(i)
            x.append(float(split_line[xind]))
            y.append(float(split_line[yind]))
            z.append(float(split_line[zind]))
            lag.append(float(split_line[lagind]))
            spect.append(float(split_line[spectind]))
            theirspect.append(timing["PMTmatch_nphe"][i-1][0])
            bubdata[evid]=[float(split_line[spectind]),float(split_line[pmttraceind])]

    i+=1
plt.figure()
plt.scatter(spect,theirspect,18)
plt.xlabel('My pulse integration',fontsize=18)
plt.ylabel('PMTmatch_nphe',fontsize=18)
plt.show

#print(bubdata)
"""

run='20170804_2'
runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
runpath = datapath+run
g = get_gain(datapath, run)
print("gain: %e"%g)
events = [evnt for evnt in os.listdir(runpath) if not os.path.isfile(os.path.join(runpath,evnt))]
Nevents = len(events)

pmtfilename = runreconpath+"PMTpulseAnalysis_%s.bin"%run
pmt = sbc.DataHandling.ReadBinary.ReadBlock(pmtfilename)

evdone = []
myphe = [0 for i in range(len(pmt["PMT_pulse_area"]))]
my_Vmax = [0 for i in range(len(pmt["PMT_pulse_area"]))]

wNIM = [0 for i in range(len(pmt["PMT_pulse_area"]))]

theirphe = [0 for i in range(len(pmt["PMT_pulse_area"]))]
ratio = [0 for i in range(len(pmt["PMT_pulse_area"]))]
j=0
for n in range(1):
    print("event %d"%n)
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath, n)
    print(len(e["PMTtraces"]["traces"]))
    for k in range(len(e["PMTtraces"]["traces"])):
        trace = np.fabs(e["PMTtraces"]["traces"][k][0]) 
        trace = trace[:-50]
        
        #if ch0 saturated, stitch in low res channel:
        if max(trace) == 128:
            trace = pi.stitchTraces(trace, np.fabs(e["PMTtraces"]["traces"][k][1]))
        dt = e["PMTtraces"]["dt"][k][0]
        baseline = np.mean(trace[:50])
        trace -= baseline
        my_Vmax[j] = max(trace)
        #integrate and convert to phe:
        [a,npeaks,totInt,pktimes] = pi.SBC_pulse_integrator_bressler(trace,dt)
        
        caught = False
        if min(e["PMTtraces"]["traces"][k][1]) < -0.5:
            caught = True

        myphe[j]=a/g
        
        if npeaks >= 1 and a/6 < 2 and not caught:
            print("event %d trace %d (had %d peaks)"%(n, k, npeaks))
            baseline_std = np.std(trace[0:50])
            tPMT = np.arange(len(trace))*dt
            integral_t0_index = np.argmax(np.diff(trace)>4)
            integral_t0 = tPMT[integral_t0_index]
            p,t = pi.get_pulse(trace,tPMT,dt, 500e-9, integral_t0,baseline_std)
            ret = scipy.integrate.trapz(p,t)*CONVERSION_TO_CHARGE
            print("ret: %f"%ret)
            print("a: %f"%a)
            
            plt.figure()
            plt.scatter(t,p, 10, c='r')
            plt.plot(tPMT,trace, c='b', lw=1)
            plt.plot([integral_t0, integral_t0], [0, max(trace)], 'g')
            plt.xlim([0, 7.5e-7])
            plt.title('%e phe'%(a/g))
            plt.plot()
            
        
        j+= 1
        gc.collect()
        

#ft = np.polyfit(myphe,ratio,0)
#print(ft)
plt.figure()
plt.scatter(myphe, my_Vmax)
plt.title('my phe vs vmax')
plt.xlabel('photoelectrons')
plt.ylabel('trace maximum (ADC)')
plt.show()    

plt.figure()
plt.hist2d(myphe, my_Vmax,150)
plt.title('my phe vs vmax')
plt.xlabel('photoelectrons')
plt.ylabel('trace maximum (ADC)')
plt.show()    


"""
run='20170628_9'
runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
runpath = datapath+run
events = [evnt for evnt in os.listdir(runpath) if not os.path.isfile(os.path.join(runpath,evnt))]
Nevents = len(events)

pmtfilename = runreconpath+"PMTpulseAnalysis_%s.bin"%run
pmt = sbc.DataHandling.ReadBinary.ReadBlock(pmtfilename)
#print(len(pmt["PMT_pulse_area"]))
print(pmt.keys())
#evid = pmt["runid"]+"-"+pmt["ev"]
#area = pmt["PMT_pulse_area"]
#bubdata[]
evdone = []
myphe = [0 for i in range(len(pmt["PMT_pulse_area"]))]
my_Vmax = [0 for i in range(len(pmt["PMT_pulse_area"]))]
theirphe = [0 for i in range(len(pmt["PMT_pulse_area"]))]
ratio = [0 for i in range(len(pmt["PMT_pulse_area"]))]
j=0
for n in range(1):
    print("event %d"%n)
    indices = [i for i,x in enumerate(pmt["ev"]) if int(x)==n]
    #print(indices)
    areas=pmt["PMT_pulse_area"][indices]
    #print(len(areas))
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,n)
    print(len(e["PMTtraces"]["traces"]))
    for k in range(len(e["PMTtraces"]["traces"])):
        trace = np.fabs(e["PMTtraces"]["traces"][k][0]) 
        #if ch0 saturated, stitch in low res channel:
        if max(trace) == 128:
            trace = pi.stitchTraces(trace,np.fabs(e["PMTtraces"]["traces"][k][1]))
        dt = e["PMTtraces"]["dt"][k][0]
        my_Vmax[k] = max(trace)/e["PMTtraces"]["v_scale"][k,0,None]
        #integrate and convert to phe:
        [a,n,totInt,pktimes] = pi.SBC_pulse_integrator_bressler(trace,dt) 
        myphe[j]=a
        theirphe[j]=-areas[k]*CONVERSION_TO_CHARGE/e["PMTtraces"]["v_scale"][k,0,None]
        ratio[j] = myphe[j]/theirphe[j]

        j+= 1
        gc.collect()
        
print(len(myphe))
print(j)
#ft = np.polyfit(myphe,ratio,0)
#print(ft)

plt.figure()
plt.scatter(myphe, my_Vmax)
plt.title('my phe vs vmax')
plt.xlabel('photoelectrons')
plt.ylabel('trace maximum (V)')
plt.show()    

plt.figure()
plt.scatter(myphe,theirphe)
plt.xlabel("My Area [elecrons]")
plt.ylabel("Their area [electrons]")
plt.xlim([min(myphe),max(myphe)])
plt.ylim([min(theirphe),max(theirphe)])
plt.show()

plt.figure()
plt.scatter(myphe,ratio)
#plt.plot([min(myphe),max(myphe)],[ft[0],ft[0]])
plt.show()

plt.figure()
plt.hist(theirphe)
plt.show()
    
"""