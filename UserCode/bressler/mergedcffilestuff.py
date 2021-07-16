#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 13:07:04 2019

@author: bressler
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import SBCcode as sbc
from dogammasmakebubbles import getRate
import runlistscatalogue as rlc
print("modules loaded")

bgruns = rlc.bgJune23and24

srcruns = rlc.cfJuly6to11

srctotbub,srctotlt,_,_,srcr,srcrerr,_=getRate(sbc.DataHandling.ReadBinary.ReadBlock("/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"),
              srcruns,5)
bgtotbub,bgtotlt,_,_,bgr,bgrerr,_=getRate(sbc.DataHandling.ReadBinary.ReadBlock("/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"),
              bgruns,5)

with open("/nashome/b/bressler/sbcoutput/bgOct10and11_merged.txt","r") as fin:
    bgdata = fin.readlines()

with open("/nashome/b/bressler/sbcoutput/YBeSept26to28_merged.txt","r") as fin:
    data = fin.readlines()


headers = data[0].split()
xind = headers.index("x")
yind = headers.index("y")
zind = headers.index("z")
at00ind = headers.index("at0_0")
at01ind = headers.index("at0_1")
pmtt0ind = headers.index("PMTt0")
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
        if float(split_line[blockedind])>0.9:
            print(split_line[blockedind])
        if (not np.isnan(float(split_line[xind]))) and float(split_line[blockedind])==0:
            goodsrcbubbles += 1
            x.append(float(split_line[xind]))
            y.append(float(split_line[yind]))
            z.append(float(split_line[zind]))
            lag.append(float(split_line[lagind]))
            spect.append(float(split_line[spectind]))

    i+=1
spect = [spect[i] for i in range(len(spect)) if not np.isnan(spect[i])]
print(spect)
bgspect = []
goodbgbubbles = 0
j=0
for line in bgdata:
    if j>0:
        split_line = line.split()
        if (not np.isnan(float(split_line[xind]))):
            goodbgbubbles += 1
            bgspect.append(float(split_line[spectind]))

    j+=1

bgspect = [bgspect[i] for i in range(len(bgspect)) if not np.isnan(bgspect[i])]

srcbubswithscintillation = len([lag[i] for i in range(len(lag)) if not np.isnan(lag[i])])
bgbubswithscintillation = len([bgspect[i] for i in range(len(bgspect)) if not np.isnan(bgspect[i])])

r2 = [x[i]**2 +y[i]**2 for i in range(len(x))]
plt.figure()
plt.grid()
for i in range(len(x)):
    if not np.isnan(lag[i]):
        plt.scatter(x[i],y[i],40,'r','s')
    else:
        plt.scatter(x[i],y[i],40)
plt.xlabel("x [cm]",fontsize=18)
plt.ylabel('y [cm]',fontsize=18)
plt.show

rjar = max(x)

plt.figure()
plt.grid()
for i in range(len(r2)):
    if not np.isnan(lag[i]):
        plt.scatter(r2[i]/rjar,z[i],40,'r','s')
    else:
        plt.scatter(r2[i]/rjar,z[i],40)
plt.xlabel("r^2/r_jar [cm]",fontsize=18)
plt.ylabel("z [cm]",fontsize=18)
plt.show

bins = [2**i+0.5 for i in range(11)]
#bins = np.arange(0.5,1+np.ceil(max(spect)))
bins = np.insert(bins,0,0.5)
bins=np.insert(bins,0,0)
bins=np.insert(bins,0,-1)
print(bins)
plt.figure()
Nsrc,_,_ = plt.hist(spect,bins,histtype='step',color='r',label="lower threshold")
Nbg,_,_=plt.hist(bgspect,bins,histtype='step',color='b',label="higher threshold")
plt.legend()
plt.show

Nsrc[0] += goodsrcbubbles - srcbubswithscintillation
Nbg[0] += goodbgbubbles - bgbubswithscintillation

phepts = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
fsrc =[Nsrc[i]/goodsrcbubbles for i in range(len(Nsrc))]
errsrch = [(0.5+np.sqrt(Nsrc[i]+0.25))/goodsrcbubbles for i in range(len(Nsrc))]
errsrcl = [(-0.5+np.sqrt(Nsrc[i]+0.25))/goodsrcbubbles for i in range(len(Nsrc))]

fbg = [Nbg[i]/goodbgbubbles for i in range(len(Nbg))]
errbgh = [(0.5+np.sqrt(Nbg[i]+0.25))/goodbgbubbles for i in range(len(Nbg))]
errbgl = [(-0.5+np.sqrt(Nbg[i]+0.25))/goodbgbubbles for i in range(len(Nbg))]



plt.figure()
plt.errorbar(phepts,fsrc,[errsrcl,errsrch],c='r',capsize=0,linestyle="",label="BiBe")
#plt.scatter(phepts,fcf,30,c='r')
plt.errorbar(phepts,fbg,[errbgl,errbgh],c='b',capsize=0,linestyle="",label="BiAl")
#plt.scatter(phepts,fbg,30,c='b')

for j in range(len(fsrc)):
    plt.plot([bins[j],bins[j]],[0,fsrc[j]],c='r')
    plt.plot([bins[j],bins[j+1]],[fsrc[j],fsrc[j]],c='r')
    plt.plot([bins[j+1],bins[j+1]],[0,fsrc[j]],c='r')
    
for j in range(len(fbg)):
    plt.plot([bins[j],bins[j]],[0,fbg[j]],c='b')
    plt.plot([bins[j],bins[j+1]],[fbg[j],fbg[j]],c='b')
    plt.plot([bins[j+1],bins[j+1]],[0,fbg[j]],c='b')
plt.xlabel("Scintillation Signal [Photoelectrons]",fontsize=18)
plt.ylabel("Fraction of Bubbles",fontsize=18)
plt.grid()
plt.ylim([0,0.3])
plt.xlim([-1,40])
plt.legend(fontsize=18)
plt.show

print(Nsrc)

binsizes = [(bins[i+1]-bins[i]) for i in range(len(Nsrc))]

def fancyerrorbar(n,h):
    if h==1:
        return 0.5 + np.sqrt(n + 0.25)
    else:
        return -0.5 + np.sqrt(n + 0.25)

rsrc=[Nsrc[i]/(srctotlt) for i in range(len(Nsrc))]
errrsrch = [fancyerrorbar(Nsrc[i],1)/(srctotlt) for i in range(len(Nsrc))]
errrsrcl = [fancyerrorbar(Nsrc[i],0)/(srctotlt) for i in range(len(Nsrc))]


rbg = [Nbg[i]/(bgtotlt) for i in range(len(Nbg))]
errrbgh = [fancyerrorbar(Nbg[i],1)/(bgtotlt) for i in range(len(Nbg))]
errrbgl = [fancyerrorbar(Nbg[i],0)/(bgtotlt) for i in range(len(Nbg))]

plt.figure()
plt.errorbar(phepts,rsrc,[errrsrcl,errrsrch],c='r',marker='^',capsize=0,linestyle="",label="YBe")
#plt.scatter(phepts,fcf,30,c='r')
plt.errorbar(phepts,rbg,[errrbgl,errrbgh],c='b',marker='v',capsize=0,linestyle="",label="background")
#plt.scatter(phepts,fbg,30,c='b')
"""
for j in range(len(fsrc)):
    plt.plot([bins[j],bins[j]],[0,rsrc[j]],c='r')
    plt.plot([bins[j],bins[j+1]],[rsrc[j],rsrc[j]],c='r')
    plt.plot([bins[j+1],bins[j+1]],[0,rsrc[j]],c='r')
    
for j in range(len(fbg)):
    plt.plot([bins[j],bins[j]],[0,rbg[j]],c='b')
    plt.plot([bins[j],bins[j+1]],[rbg[j],rbg[j]],c='b')
    plt.plot([bins[j+1],bins[j+1]],[0,rbg[j]],c='b')
"""
plt.xlabel("Scintillation Signal [Photoelectrons]",fontsize=18)
plt.ylabel("Rate [Hz/bin]",fontsize=18)
plt.grid()
plt.ylim([0,0.05])
plt.xlim([0,1e3])
plt.xscale('symlog',linthreshx=1)
plt.yscale('symlog',linthreshy=1e-7)
plt.legend(fontsize=18)
plt.show


cfcdf = np.cumsum(fsrc)
bgcdf = np.cumsum(fbg)

plt.figure()
plt.plot(phepts,cfcdf)
plt.plot(phepts,bgcdf)
plt.xlabel("Photoelectrons")
plt.ylabel("CDF")
plt.show

print(1.858*np.sqrt((goodbgbubbles+goodsrcbubbles)/(goodbgbubbles*goodsrcbubbles)))

[D,p] = stats.ks_2samp(fsrc,fbg)
print(D)
print(p)

print(srcbubswithscintillation/goodsrcbubbles)
print(bgbubswithscintillation/goodbgbubbles)