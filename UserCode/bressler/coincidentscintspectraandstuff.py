#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:56:08 2019

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import matplotlib.pyplot as plt
import scipy
from runlistscatalogue import *

bgruns = bgOct10and11
cfruns = BiAlOct6to9

allruns = [bgruns,cfruns]

bgspect = []
cfspect = []

bgzs = []
cfzs = []

bglags = []
cflags = []
j=0
for runs in allruns:
    ch0files = []
    ch1files = []
    lines = []
    datas = []
    for run in runs:
        ch0files.append(open("/nashome/b/bressler/sbcoutput/%s_PMTmatching_ch0.txt"%run,"r"))
        #ch1files.append(open("/nashome/b/bressler/output/%s_PMTmatching_ch1.txt"%run,"r"))
    
    for f in ch0files:
    
        for i,l in enumerate(f):
            if i>=1:
                d=l.split()
                run_num = d[0]
                event_num = d[1]
                ind = d[2]
                lag = d[3]
                pmtt0 = d[4]
                at0 = d[5]
                phe = d[6]
                z = d[7]
                if j == 0:
                    bgspect.append(float(phe))
                    bgzs.append(float(z))
                    bglags.append(float(lag))
                elif j == 1:
                    cfspect.append(float(phe))
                    cfzs.append(float(z))
                    cflags.append(float(lag))

        f.close()
    j += 1
plt.figure()
plt.grid()
plt.scatter(bgzs,bgspect,20,'r',label="background")
plt.scatter(cfzs,cfspect,20,'g',label="californium")
plt.legend(fontsize=18)
plt.xlabel("z position (cm)",fontsize=18)
plt.ylabel("scintillation signal (phe)",fontsize=18)
plt.show

plt.figure()
plt.grid()
plt.scatter(bglags,bgspect,20,'r',label="background")
plt.scatter(cflags,cfspect,20,'g',label="californium")
plt.legend(fontsize=18)
plt.xlabel("lag (microseconds)",fontsize = 18)
plt.ylabel("scintillation signal (phe)",fontsize = 18)
plt.show


plt.figure()
plt.grid()
print(int(np.ceil(max(bgspect))))
bins = np.arange(0,int(max(bgspect)+1),1)
bgN,_,_ = plt.hist(bgspect,bins,histtype = "step")
print(bins)
bgphepts = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
plt.errorbar(bgphepts,bgN,np.sqrt(bgN))
plt.show


plt.figure()
plt.grid()
bins = np.arange(0,int(max(cfspect)+1),1)
print(bins)
cfN,_,_ = plt.hist(cfspect,bins,histtype = "step")
cfphepts = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
plt.errorbar(cfphepts,cfN,np.sqrt(cfN))
plt.show

plt.figure()
plt.grid()
plt.hist(cfspect,bins,histtype="step",normed=True)
plt.hist(bgspect,bins,histtype="step",normed=True)
plt.show

plt.figure()
plt.grid()
plt.scatter(bgphepts,bgN/len(bgspect),40,c='r',label="background")
plt.scatter(cfphepts,cfN/len(cfspect),40,c='g',label="californium")
plt.errorbar(bgphepts,bgN/len(bgspect),np.sqrt(bgN)/len(bgspect),capsize=0,ls="")
plt.errorbar(cfphepts,cfN/len(cfspect),np.sqrt(cfN)/len(cfspect),capsize=0,ls="")
plt.legend(fontsize=18)
plt.show