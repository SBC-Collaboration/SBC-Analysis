#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:31:59 2020

@author: bressler
"""

import SBCcode as sbc
import matplotlib.pyplot as plt
import runlistscatalogue as rlc
import numpy as np
from os.path import isfile,join
from os import listdir
import gc

bial = rlc.BiAlOct6to9
cf = rlc.cfJuly6to11
bg = rlc.bgOct10and11
bg2 = rlc.bgOct2and3
bg3 = rlc.bgCombinedMultiTemp
types = [cf]
runs = []
lt = 0
for typE in types:
    for run in typE:
        runs.append(run)
evids = []
phe = []
for run in runs:
    print(run)
    with open("/nashome/b/bressler/sbcoutput/%s_muonCoincidences.txt"%run,"r") as coincfile:
        data = coincfile.readlines()
    if len(data)>1:
        for i in range(1,len(data)):
            evids.append(data[i].split()[0]+"-"+data[i].split()[1])
            phe.append(float(data[i].split()[2].rstrip()))
            if phe[-1]<1:
                print(evids[-1])
    runpath = '/bluearc/storage/SBC-17-data/'+run
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    
    for event in events:
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
        lt += e["event"]["livetime"]
        gc.collect()
print(evids)
print(phe)
print(len(phe))
print(lt)
print("Rate: %f Hz"%(len(phe)/lt))
bins = [(2**i)+0.5 for i in range(12)]
#bins = np.arange(1,int(1+np.ceil(max(spect))))
bins = np.insert(bins,0,0.5)
bins=np.insert(bins,0,-0.5)
bins=np.insert(bins,0,-1.5)
binc=[(bins[i+1]+bins[i])/2 for i in range(len(bins)-1)]
plt.figure()
N,_,_=plt.hist(phe,bins,histtype='step',linewidth=4)
plt.xscale('symlog')
plt.grid()
plt.xlabel('phe')
plt.ylabel('count')
plt.show()
plt.figure()
plt.errorbar(binc,N,np.sqrt(N),ds='steps-mid')
plt.xscale('symlog')
plt.xlabel('phe')
plt.show