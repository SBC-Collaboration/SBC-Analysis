#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 14:12:39 2019

@author: bressler
"""

import scipy
import numpy as np
import matplotlib.pyplot as plt
from coincidentbubblescintillation import zdependence


bgruns = ["20170711_17","20170711_18","20170712_0","20170712_1"]

cfruns = ["20170710_8","20170710_9"]

"""
        
        "20170707_6","20170707_7","20170707_8","20170707_9","20170707_10","20170708_0",
        "20170708_1","20170708_3","20170708_4","20170708_5","20170708_6",
        "20170708_7","20170708_8","20170708_9","20170709_0","20170709_1","20170709_2",
        "20170709_3","20170709_4","20170709_6","20170709_7","20170709_8",
        ,"20170711_0","20170711_14","20170711_15","20170711_16"]

        files made:
            "20170710_0","20170710_1","20170710_2","20170710_3","20170710_4","20170710_5",
            "20170710_6","20170710_7",

bad AcousticAnalysis_ .bin files:
    "20170708_2",
"""

#pmtnobubdiffs,pmtdiffs,dubbubdiffs = trig_difference(cfruns)
goodz,diffs,coincspec,Ncoinc,ntotcoinc,totevents,totbub = zdependence(cfruns)

#bgpmtnobubdiffs,bgpmtdiffs,bgdubbubdiffs = trig_difference(bgruns)
#bggoodz,bgdiffs,bgcoincspec,bgNcoinc,bgntotcoinc,bgtotevents,bgtotbub = zdependence(bgruns)

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
#bgvals,_,_ = plt.hist(bgcoincspec,bins,histtype='step')
plt.xlabel("spectrum of PMT pulse areas within 500 microseconds before acoustic t_0 (photoelectrons)")
plt.show

plt.figure()
plt.bar(bins[:(len(vals))],[v/totbub for v in vals],1,color='r',linewidth=0,label="californium")
#plt.bar(bins[:len(vals)],[v/bgtotbub for v in bgvals],0.7,color='b',linewidth = 0,label="background")
plt.xlabel("spectrum of PMT pulse areas within 500 microseconds before acoustic t_0 (photoelectrons)",fontsize=25)
plt.ylabel("probability of a scintillation pulse of this area per bubble",fontsize=25)
plt.legend(fontsize=18)
plt.show
