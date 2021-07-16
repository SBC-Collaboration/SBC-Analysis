#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:56:08 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
import runlistscatalogue as rlc


runs = ["20170707_6","20170707_7","20170707_8","20170707_9","20170707_10","20170708_0",
            "20170708_1","20170708_3","20170708_4","20170708_5","20170708_6",
            "20170708_7","20170708_8","20170708_9","20170709_0","20170709_1","20170709_2",
            "20170709_3","20170709_4"]
ch0files = []
ch1files = []
lines = []
ch0zs = []
ch0lags = []
ch1zs = []
ch1lags = []
spect = []
allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
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
            spect.append(float(phe))
            ch0zs.append(float(z))
            ch0lags.append(float(lag))


    f.close()
"""   
for f in ch1files:

    for i,l in enumerate(f):
        if i>=1:
            d=l.split()
            run_num = d[0]
            event_num = d[1]
            ind = d[2]
            lag = d[3]
            pmtt0 = d[4]
            at0 = d[5]
            z = d[7]
            ch1zs.append(float(z))
            ch1lags.append(float(lag))


    f.close()
"""
def fit_fun(x, A, k, ph, m, b):
    return m*x + b + A*np.cos(k*x + ph)

def gaussian(x,mu,sigma,amplitude):
    return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )

def lin(x,m,b):
    return m*x + b

zbinedges = np.arange(-3,0.5,0.5)
zpts = [(zbinedges[i]+zbinedges[i+1])/2 for i in range(len(zbinedges)-1)]
zbinlags = [[] for i in range(len(zbinedges)-1)]
print(zbinlags)
for j in range(len(ch0lags)):
    for i in range(len(zbinedges)-1):
        if ch0zs[j] >= zbinedges[i] and ch0zs[j] <zbinedges[i+1]:
            zbinlags[i].append(ch0lags[j])
"""
params, params_cov = scipy.optimize.curve_fit(fit_fun, ch1zs, ch1lags)

fitx = np.arange(-3,0,0.1)
ft = [fit_fun(fitx[i], params[0], params[1], params[2],
              params[3], params[4]) for i in range(len(fitx))]

ss = (-1/params[3])*10000
"""

avgs = []

gausfitx = np.arange(-400,0,0.1)

for j in range(len(zpts)):
    Nbins=int(np.floor(np.sqrt(len(zbinlags[j]))))
    print(Nbins)
    if Nbins == 0:
        Nbins = 1
    if len(zbinlags[j])>0:
        dat,bins = np.histogram(zbinlags[j],Nbins)
        bincenters = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        gparams,gparams_cov = scipy.optimize.curve_fit(gaussian,bincenters,dat,
                                                           p0=(-200,50,45))
        avgs.append(gparams[0])
        gausfit = [gaussian(gausfitx[i],gparams[0],gparams[1],gparams[2]) for i in range(len(gausfitx))]

    
lfitx = np.arange(-3,0,0.1)
print(zpts)
print(avgs)
lparams, lparams_cov = scipy.optimize.curve_fit(lin,zpts,avgs)
ss = (-1/lparams[0])*10000
plt.figure()
plt.grid()
plt.scatter(ch0zs,ch0lags)
plt.plot(lfitx,lparams[0]*lfitx + lparams[1], color="g",linewidth=4,
         label = "fit line: dt = %fz +%f"%(lparams[0],lparams[1]))
#plt.plot(fitx,params[3]*fitx + params[4],color = 'g',linewidth=4,
#         label = "Linear part of fit, dt = %fz+%f"%(params[3],params[4]))
#plt.plot(fitx,ft,color='r',linewidth=5,
#         label = "Fit: dt = %fcos(%fz+%f) + %fz +%f"%(params[0],
#                                 params[1],params[2],params[3],params[4]))
plt.scatter(zpts,avgs,color='c',s=100)
plt.xlabel("z position [cm]",fontsize=20)
plt.ylabel("Lag time (t_PMT - t_acoustic) [microseconds]",fontsize=20)
plt.title("Channel 0, sound speed=%s m/s"%str(ss),fontsize=20)
plt.legend(fontsize=18)
plt.show

plt.figure()
plt.hist(spect,int(np.floor(np.sqrt(len(spect)))))
plt.show
