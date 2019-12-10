#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 16:13:14 2019

@author: bressler
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

with open("/nashome/b/bressler/sbcoutput/allYBe_merged.txt",'r') as fin:
    data = fin.readlines()
    
headers = data[0].split()
xind = headers.index("x")
yind = headers.index("y")
zind = headers.index("z")
ap1ind = headers.index("AP1")
ap2ind = headers.index("AP2")
at00ind = headers.index("at0_0")
at01ind = headers.index("at0_1")
pmtt0ind = headers.index("PMTt0")
lagind = headers.index("lag")
spectind = headers.index("PMTphe")
blockedind = headers.index('isBlocked')
runind = headers.index("run")
evind = headers.index("event")
allAP2 = []
allAP1 = []
AP1 = []
AP2=[]
phe = []
x=[]
y=[]
z=[]
eventtags = []

bl = []

i=0
for line in data:
    if i>0:
        split_line = line.split()
        allAP1.append(float(split_line[ap1ind]))
        allAP2.append(float(split_line[ap2ind]))
        if not np.isnan(float(split_line[spectind])) and not np.isnan(float(split_line[xind])):
            AP1.append(float(split_line[ap1ind]))
            AP2.append(float(split_line[ap2ind]))
            phe.append(float(split_line[spectind]))
            x.append(float(split_line[xind]))
            y.append(float(split_line[yind]))
            z.append(float(split_line[zind]))
        eventtags.append(str(split_line[runind]+"-"+split_line[evind]))
        bl.append(float(split_line[blockedind]))
    i+=1
    
print(x)
    
for p in AP1:
    if p <= 0:
        ind = AP1.index(p)
        AP1.remove(p)
        AP2.remove(AP2[ind])
        phe.remove(phe[ind])
        x.remove(x[ind])
        y.remove(y[ind])
        z.remove(z[ind])
        
for p in AP2:
    if p <= 0:
        ind = AP2.index(p)
        AP2.remove(p)
        AP1.remove(AP1[ind])
        phe.remove(phe[ind])
        x.remove(x[ind])
        y.remove(y[ind])
        z.remove(z[ind])
        
r2 = [x[i]**2 + y[i]**2 for i in range(len(x))]
        
allAP1 = [ap for ap in allAP1 if ap>0]
allAP2 = [ap for ap in allAP2 if ap>0]

logAllAP1 = [np.log(x) for x in allAP1]
logAllAP2 = [np.log(x) for x in allAP2]
print(min(AP1))
logAP1 = [np.log(x) for x in AP1]
logAP2 = [np.log(x) for x in AP2]

plt.figure()
plt.hist(logAllAP1,int(np.floor(np.sqrt(len(logAllAP1)))),histtype='step')
plt.xlabel("uncorrected ln(AP1)")
plt.ylabel("count")
plt.show

plt.figure()
plt.scatter(logAP1,phe)
plt.xlabel("uncorrected ln(AP1)",fontsize=18)
plt.ylabel("Scintillation signal [photoelectrons]",fontsize=18)
plt.show

plt.figure()
plt.scatter(r2,logAP1)
plt.xlabel("r^2 (cm^2)")
plt.ylabel("uncorrected ln(AP1)")
plt.show

plt.figure()
plt.scatter(z,logAP2,c='r',label="Piezo 2")
plt.scatter(z,logAP1,c='k',label="Piezo 1")
plt.xlabel("z [cm]",fontsize=18)
plt.ylabel("Uncorrected ln(AP)",fontsize=18)
plt.legend(fontsize=18)
plt.show


def gaussian(x,mu,sigma,amplitude):
    return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )

def lin(x,m,b):
    return m*x + b

def exp(x,lam,A):
    return A*np.exp(lam*x)


zbinedges = np.linspace(-2.9,0,6)
zpts = [(zbinedges[i]+zbinedges[i+1])/2 for i in range(len(zbinedges)-1)]
zbinAP1s = [[] for i in range(len(zbinedges)-1)]
for j in range(len(logAP1)):
    for i in range(len(zbinedges)-1):
        if z[j] >= zbinedges[i] and z[j] < zbinedges[i+1] and not np.isinf(logAP1[j]):
            zbinAP1s[i].append(logAP1[j])


avgs1 = []

gausfitx = np.arange(-400,0,0.1)

for j in range(len(zpts)):
    Nbins=int(np.floor(np.sqrt(len(zbinAP1s[j]))))
    if Nbins == 0:
        Nbins = 1
    if len(zbinAP1s[j])>0:
        med = np.median(zbinAP1s[j])
        dat,bins = np.histogram(zbinAP1s[j],Nbins)
        bincenters = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        try:
            gparams,gparams_cov = curve_fit(gaussian,bincenters,dat,
                                                               p0=(-11,3,15))
            print("sigma = "+str((gparams[0]-med)/gparams[1]))
            if (gparams[0]-med)/gparams[1] > 0.5:
                avgs1.append(med)
            else:
                avgs1.append(gparams[0])
            gausfit = [gaussian(gausfitx[i],gparams[0],gparams[1],gparams[2]) for i in range(len(gausfitx))]
        except:
            print("couldn't get gaussian fit piezo 1")
            avgs1.append(med)
        
lfitx = np.arange(-3,0,0.1)
lparams1, lparams_cov1 = curve_fit(lin,zpts,avgs1)

zbinAP2s = [[] for i in range(len(zbinedges)-1)]
for j in range(len(logAP2)):
    for i in range(len(zbinedges)-1):
        if z[j] >= zbinedges[i] and z[j] < zbinedges[i+1] and not np.isinf(logAP2[j]) and logAP2[j]>-12:
            zbinAP2s[i].append(logAP2[j])


avgs2 = []

gausfitx = np.arange(-20,0,0.1)

for j in range(len(zpts)):
    Nbins=int(np.floor(np.sqrt(len(zbinAP2s[j]))))
    if Nbins == 0:
        Nbins = 1
    if len(zbinAP2s[j])>0:
        med = np.median(zbinAP2s[j])
        dat,bins = np.histogram(zbinAP2s[j],Nbins)
        bincenters = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
        try:
            gparams,gparams_cov = curve_fit(gaussian,bincenters,dat,
                                                               p0=(-5,4,15))
            if (gparams[0]-med)/gparams[1] > 0.5:
                avgs2.append(med)
            else:
                avgs2.append(gparams[0])
            gausfit = [gaussian(gausfitx[i],gparams[0],gparams[1],
                                gparams[2]) for i in range(len(gausfitx))]

        except:
            print("couldn't get gaussian fit piezo 2")

            avgs2.append(med)
        
lfitx = np.arange(-3,0,0.1)
lparams2, lparams_cov2 = curve_fit(lin,zpts,avgs2)

zcorrAP1 = np.zeros(len(logAP1))
for i in range(len(logAP1)):
    if not np.isinf(logAP1[i]):
        zcorrAP1[i] = logAP1[i]-(lparams1[0]*z[i]+lparams1[1])
        
zcorrAP2 = np.zeros(len(logAP2))
for i in range(len(logAP2)):
    if not np.isinf(logAP2[i]):
        zcorrAP2[i] = logAP2[i]-(lparams2[0]*z[i]+lparams2[1])

plt.figure()
plt.grid()
plt.scatter(z,logAP1)
plt.scatter(z,logAP2)
plt.plot(lfitx, lparams1[1]+(lparams1[0]*lfitx), color="g",linewidth=4,
         label = "fit line piezo 1")
plt.plot(lfitx,lparams2[1]+(lparams2[0]*lfitx),color="r",linewidth=4,
         label="fit line piezo 2")
#plt.plot(fitx,params[3]*fitx + params[4],color = 'g',linewidth=4,
#         label = "Linear part of fit, dt = %fz+%f"%(params[3],params[4]))
#plt.plot(fitx,ft,color='r',linewidth=5,
#         label = "Fit: dt = %fcos(%fz+%f) + %fz +%f"%(params[0],
#                                 params[1],params[2],params[3],params[4]))
plt.scatter(zpts,avgs1,color='c',s=100)
plt.scatter(zpts,avgs2,color='m',s=100)
plt.xlabel("z position [cm]",fontsize=20)
plt.ylabel("ln(AP1)",fontsize=20)
plt.legend(fontsize=18)
plt.show



    
plt.figure()
plt.scatter(zcorrAP1,phe)
plt.xlabel("corrected ln(AP1)",fontsize=18)
plt.ylabel("scintillation signal (pe)",fontsize=18)
plt.yscale('log')
plt.grid()
plt.show

plt.figure()
plt.scatter(z,zcorrAP1)
plt.xlabel("z (cm)",fontsize=18)
plt.ylabel("corrected ln(AP1)",fontsize=18)
plt.show

plt.figure()
plt.hist(zcorrAP1,int(np.floor(np.sqrt(len(zcorrAP1)))),histtype='step',label="piezo 1")
plt.hist(zcorrAP2,int(np.floor(np.sqrt(len(zcorrAP2)))),histtype='step',label="piezo 2")

plt.yscale('log')
plt.xlabel("corrected ln(AP)",fontsize=18)
plt.legend(fontsize=18)
plt.show

plt.figure()
plt.scatter(zcorrAP1,zcorrAP2,50)
plt.show

combinedAP = [(zcorrAP1[i]+zcorrAP2[i])/2 for i in range(len(zcorrAP1))]
plt.figure()
plt.hist(combinedAP,int(np.floor(np.sqrt(len(combinedAP)))),histtype='step')
plt.xlabel('combined ln(AP)',fontsize=18)
plt.ylabel('counts per bin')
plt.show

for i in range(len(combinedAP)):
    if combinedAP[i]<-2:
        print(eventtags[i])