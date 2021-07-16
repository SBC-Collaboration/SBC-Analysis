#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:01:02 2019

@author: bressler
"""

import matplotlib.pyplot as plt
plt.style.use('default')
from matplotlib.lines import Line2D
import numpy as np
import scipy

import gc
gc.collect()

# -50C data:
Q_bg50 = [2.00, 2.16]
r_bg50 = [2.81, 2.12]
r_err_bg50 = [0.38, 0.32]

Q_cf50 = [1.91,2.08]
r_cf50 = [10.01, 7.34]
r_err_cf50 = [3.64, 0.36]

#-45C data

Q_bg45a = [1.14] # july 11-12
r_bg45a = [2.53]
r_err_bg45a = [0.26]

Q_bg45b = [1.14, 1.33, 1.44, 1.56, 1.85] # sept 26
r_bg45b = [3.54, 3.24, 1.94, 4.95, 2.71]
r_err_bg45b = [1.36, 1.25, 1.45, 1.49, 1.13]

Q_bg45c = [1.14, 1.33, 1.44, 1.56] #sept 27
r_bg45c = [4.81, 3.76, 2.75, 2.02]
r_err_bg45c = [1.17, 1.00, 0.88, 0.73]

Q_YBeO45 = [1.14, 1.33, 1.44, 1.56]
r_YBeO45 = [15.84, 12.27, 11.22, 11.33]
r_err_YBeO45 = [2.24, 1.74, 1.65, 1.57]

Q_cf45 = [1.14]
r_cf45 = [11.42]
r_err_cf45 = [0.4] 

Q_BiAl45 = [1.14]
r_BiAl45 = [4.01]
r_err_BiAl45 = [0.56]

#Q_Bi45 = [1.14]
#r_Bi45 = [8.22]
#r_err_Bi45 = [4.93]

#Q_co45a = [1.14]
#r_co45a = [4.35]
#r_err_co45a = [3.24]

Q_co45b = [1.14, 1.33, 1.43, 1.56]
r_co45b = [6.78, 4.57, 2.91, 6.38]
r_err_co45b = [3.49, 2.74, 2.18, 3.28]


Q_BiBe45 = [1.14]
r_BiBe45 = [3.14]
r_err_BiBe45 = [0.48]

Q_yal45 = [1.14, 1.33, 1.44, 1.56]
r_yal45 = [4.56, 4.63, 3.96, 3.11]
r_err_yal45 = [1.11, 1.09, 0.99, 0.83]

#-43C data:

Q_bg43a = [0.90] # august 3 and 4
r_bg43a = [4.80]
r_err_bg43a = [0.38]

Q_bg43b = [0.90, 0.97, 1.04, 1.12] # October 2 and 3
r_bg43b = [3.27, 4.48, 3.25, 4.01]
r_err_bg43b = [0.54, 1.06, 0.81, 0.97]

Q_ybe43 = [0.90, 0.97, 1.04, 1.12]
r_YBeO43 = [15.38, 18.72, 16.86, 12.19]
r_err_YBeO43 = [3.36, 3.12, 2.98, 2.44]

Q_cf43 = [0.90]
r_cf43 = [12.08]
r_err_cf43 = [1.31] 

Q_bibe43a = [0.90] # July 13
r_BiBe43a = [3.97]
r_err_BiBe43a = [0.43]

Q_bibe43b = [0.97, 1.04, 1.12] # Sept 30
r_BiBe43b = [2.15, 6.21, 5.37]
r_err_BiBe43b = [0.99, 1.72, 1.62]

Q_co43 = [0.90, 0.97, 1.04, 1.12]
r_co43 = [6.14, 4.05, 1.39, 1.40]
r_err_co43 = [2.57, 2.43, 1.39, 1.40]


# -41C
Q_bibe41 = [0.72]
r_bibe41 = [4.80]
r_err_bibe41 = [0.33]

Q_bial41 = [0.72]
r_bial41 = [5.14]
r_err_bial41 = [0.39]

Q_yal41 = [0.71, 0.76, 0.81, 0.87]
r_yal41 = [3.39, 7.60, 8.36, 5.12]
r_err_yal41 = [1.74, 4.56, 4.30, 3.07]

#-38C
Q_bg38 = [0.51, 0.53, 0.57, 0.60]
r_bg38 = [12.21, 10.93, 8.96, 11.38]
r_err_bg38 = [4.40, 5.65, 5.38, 5.85]

Q_yal38 = [0.5, 0.57, 0.6]
r_yal38 = [19.38, 7.20, 9.89]
r_err_yal38 = [11.63, 7.20, 9.89]

Q_co38 = [0.56]
r_co38 = [5.93]
r_err_co38 = [4.43]

# all backgrounds, to fit:
Qbglists = [Q_bg45a, Q_bg45b, Q_bg45c, Q_bg43a, Q_bg43b, Q_bg38]
r_bglists = [r_bg45a, r_bg45b, r_bg45c, r_bg43a, r_bg43b, r_bg38]
r_err_bglists = [r_err_bg45a, r_err_bg45b, r_err_bg45c, r_err_bg43a, r_err_bg43b, r_err_bg38]

# all YBe to fit:
Qybelists = [Q_YBeO45, Q_ybe43]
r_ybelists = [r_YBeO45, r_YBeO43]
r_err_ybelists = [r_err_YBeO45,r_err_YBeO43]

Q_allbg = []
r_allbg = []
rerr_allbg = []

for i in range(len(Qbglists)):
    q = Qbglists[i]
    r = r_bglists[i]
    rerr = r_err_bglists[i]
    for j in range(len(q)):
        Q_allbg.append(q[j])
        r_allbg.append(r[j])
        rerr_allbg.append(rerr[j])
        
Q_allybe = []
r_allybe = []
rerr_allybe = []

for i in range(len(Qybelists)):
    q = Qybelists[i]
    r = r_ybelists[i]
    rerr = r_err_ybelists[i]
    for j in range(len(q)):
        Q_allybe.append(q[j])
        r_allybe.append(r[j])
        rerr_allybe.append(rerr[j])
        
linex = np.arange(0.5,2.5,0.1)

#fitfunc = lambda p, x: p[0] + p[1] * np.array(x)
flatfitfunc = lambda p, x: p[0]

#errfunc = lambda p, x, y, err: (y - fitfunc(p,x))/err
flaterrfunc = lambda p, x, y, err: (y - flatfitfunc(p,x))/err

pinit = [1]
exppinit = [0,0.1,-1]
#out = scipy.optimize.leastsq(errfunc, pinit,
#                             args = (Q_allybe, r_allybe, rerr_allybe), full_output=1)
flatout = scipy.optimize.leastsq(flaterrfunc, exppinit,
                             args = (Q_allbg, r_allbg, rerr_allbg), full_output=1)
flatpfinal = flatout[0]
print(flatpfinal)
#pfinal = out[0]
#cov = out[1]
#print(pfinal)
#print(cov)

def line(x,m,b):
    return m*x + b

#expy = exppfinal[0] + exppfinal[1]*np.exp(exppfinal[2]*linex)
flaty = [flatpfinal[0] for i in range(len(linex))]
print(linex)
qsubtract = 0.71
print("%f keV - > %f mHz"%(qsubtract,flatpfinal[0]))
bgresid = np.zeros([len(r_allbg),1])
for i in range(len(bgresid)):
    bgresid[i] = r_allbg[i] - flatfitfunc(flatpfinal,Q_allbg[i])
bgresidstd = np.std(bgresid)
print("bg residual standard deviation: %f mHz"%bgresidstd)
flaty1 = flaty + 1*bgresidstd
flaty2 = flaty - 1*bgresidstd

#lineyscipy = pfinal[0] + pfinal[1]*linex
#y1 = line(linex,pfinal[1]+cov[1,1]**0.5, pfinal[0]-cov[0,0]**0.5)
#y2 = line(linex,pfinal[1]-cov[1,1]**0.5, pfinal[0]+cov[0,0]**0.5)

#ybeexpfitfunc = lambda p,x: exppfinal[0] + p[0] * (np.exp(p[1]*np.array(x)))
#ybeexperrfunc = lambda p, x, y, err: (y - ybeexpfitfunc(p,x))/err

#ybeexppinit = [10,-1]

#ybeexpout = scipy.optimize.leastsq(ybeexperrfunc,ybeexppinit,
#                             args = (Q_allybe, r_allybe, rerr_allybe), full_output=1)
#ybeexppfinal = ybeexpout[0]
#ybeexpy =exppfinal[0] + ybeexppfinal[0]*np.exp(ybeexppfinal[1]*linex)

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(linex,lineyscipy, linewidth = 3, label='Linear fit to Y-Be')
ax.plot(linex, flaty, linewidth=3, label='constant fit to background')
#ax.plot(linex, ybeexpy, linewidth=3, label='exponential fit to Y-Be')

ax.plot(linex,flaty1,'k--')
ax.plot(linex,flaty2,'k--')
ax.fill_between(linex,flaty1,flaty2,facecolor='gray',alpha=0.1)
elw = 2

ax.errorbar(Q_bg38,r_bg38,r_err_bg38,linestyle='none',color='r',capsize=0,elinewidth=elw,
             marker='o',label = "background -38C")
ax.errorbar(Q_bg43a,r_bg43a,r_err_bg43a,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='o',label="background -43C, 25 psia only")
ax.errorbar(Q_bg43b,r_bg43b,r_err_bg43b,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='o',label="background -43C, scan")
ax.errorbar(Q_bg45a,r_bg45a,r_err_bg45a,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='o',label="background -45C, 25 psia only")
ax.errorbar(Q_bg45b,r_bg45b,r_err_bg45b,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='o',label="background -45C, scan")
ax.errorbar(Q_bg45c,r_bg45c,r_err_bg45c,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='o',label="background -45C, scan")

ax.errorbar(Q_bibe41,r_bibe41,r_err_bibe41,linestyle='none',color='m',capsize=0,elinewidth=elw,
             marker='s',label = "Bi-Be -41C")
ax.errorbar(Q_BiBe45,r_BiBe45,r_err_BiBe45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='s',label = "Bi-Be -45C")
ax.errorbar(Q_bibe43a,r_BiBe43a,r_err_BiBe43a,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='s',label = "Bi-Be -43C, 25 psia only")
ax.errorbar(Q_bibe43b,r_BiBe43b,r_err_BiBe43b,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='s',label = "Bi-Be -43C, scan")

ax.errorbar(Q_YBeO45,r_YBeO45,r_err_YBeO45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='D',label="Y-Be -45C")
ax.errorbar(Q_ybe43,r_YBeO43,r_err_YBeO43,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='D',label="Y-Be -43C")

ax.errorbar(Q_cf43,r_cf43,r_err_cf43,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='P',label="Cf -43C")
ax.errorbar(Q_cf45,r_cf45,r_err_cf45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='P',label="Cf -45C")
ax.errorbar(Q_cf50,r_cf50,r_err_cf50,linestyle='none',color='k',capsize=0,elinewidth=elw,
             marker='P',label="Cf -50C")

ax.errorbar(Q_bial41,r_bial41,r_err_bial41,linestyle='none',color='m',capsize=0,elinewidth=elw,
             marker='<',markersize=10,label = "Bi-Al -41C")
ax.errorbar(Q_BiAl45,r_BiAl45,r_err_BiAl45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='<',markersize=10,label = "Bi-Al -45C")
#ax.errorbar(Q_Bi45,r_Bi45,r_err_Bi45,linestyle='none',color='b',capsize=0,elinewidth=elw,
#             marker='<',markersize=10,label = "Bi-207 -45C")

ax.errorbar(Q_yal45,r_yal45,r_err_yal45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='^',markersize=10,label = "Y-Al -45C")
ax.errorbar(Q_yal41,r_yal41,r_err_yal41,linestyle='none',color='m',capsize=0,elinewidth=elw,
             marker='^',markersize=10,label = "Y-Al -41C")
ax.errorbar(Q_yal38,r_yal38,r_err_yal38,linestyle='none',color='r',capsize=0,elinewidth=elw,
             marker='^',markersize=10,label = "Y-Al -38C")

#ax.errorbar(Q_co45a,r_co45a,r_err_co45a,linestyle='none',color='b',capsize=0,elinewidth=elw,
#             marker='v',markersize=10,label = "Co-57 -45C, 25 psia only (July)")
ax.errorbar(Q_co45b,r_co45b,r_err_co45b,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='v',markersize=10,label = "Co-57 -45C")
ax.errorbar(Q_co43,r_co43,r_err_co43,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='v',markersize=10,label = "Co-57 -43C")
ax.errorbar(Q_co38,r_co38,r_err_co38,linestyle='none',color='r',capsize=0,elinewidth=elw,
             marker='v',markersize=10,label = "Co-57 -38C")

ax.set_xlabel('Seitz threshold [keV]',fontsize=20)
ax.set_ylabel('Bubble Rate [mHz]',fontsize=20)
ax.grid()

custom_lines = [Line2D([0],[0],color='k',linewidth=5),
                Line2D([0],[0],color='b',linewidth=5),
                Line2D([0],[0],color='g',linewidth=5),
                Line2D([0],[0],color='m',linewidth=5),
                Line2D([0],[0],color='r',linewidth=5)]
custom_markers = [Line2D([0],[0],color='c',marker='o',label='background'),
                  Line2D([0],[0],color='c',marker='s',label='Bi-Be'),
                  Line2D([0],[0],color='c',marker='D',label='Y-Be'),
                  Line2D([0],[0],color='c',marker='P',label='Cf-252'),
                  Line2D([0],[0],color='c',marker='v',markersize=10,label='Co-57'),
                  Line2D([0],[0],color='c',marker='^',markersize=10,label='Y-Al'),
                  Line2D([0],[0],color='c',marker='<',markersize=10,label='Bi-Al')]


leg1 = ax.legend(custom_lines, ['-50C','-45C','-43C','-41C','-38C'], loc='lower left')
leg2 = ax.legend(handles=custom_markers, loc='upper right')
ax.add_artist(leg1)
#ax.grid(which='both')
#plt.xticks(np.arange(0.4,2.2,0.1))
plt.show

fig = plt.figure()

# Background
ax1 = fig.add_subplot(221)
ax1.plot(linex, flaty, linewidth=3, label='constant fit to background')
ax1.plot(linex,flaty1,'k--',label=r'$1\sigma$ region')
ax1.plot(linex,flaty2,'k--')
ax1.fill_between(linex,flaty1,flaty2,facecolor='gray',alpha=0.1)
ax1.errorbar(Q_bg38,r_bg38,r_err_bg38,linestyle='none',color='r',capsize=0,elinewidth=elw,
             marker='o',label = "-38C")
ax1.errorbar(Q_bg43a,r_bg43a,r_err_bg43a,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='o',fillstyle='none',label="-43C, 25 psia only")
ax1.errorbar(Q_bg43b,r_bg43b,r_err_bg43b,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='o',label="-43C, scan")
ax1.errorbar(Q_bg45a,r_bg45a,r_err_bg45a,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='o',fillstyle='none',label="-45C, 25 psia only")
ax1.errorbar(Q_bg45b,r_bg45b,r_err_bg45b,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='o',label="-45C, scan")
ax1.errorbar(Q_bg45c,r_bg45c,r_err_bg45c,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='o',label="-45C, scan")
ax1.set_xlabel("Qseitz [keV]")
ax1.set_ylabel("Rate [mHz]")
ax1.set_title("Backgrounds",loc='right')
ax1.legend()
ax1.grid()

# Cf-252
ax2 = fig.add_subplot(222)
ax2.plot(linex, flaty, linewidth=3, label='constant fit to background')
ax2.plot(linex,flaty1,'k--',label=r'$1\sigma$ region')
ax2.plot(linex,flaty2,'k--')
ax2.fill_between(linex,flaty1,flaty2,facecolor='gray',alpha=0.1)
ax2.errorbar(Q_cf43,r_cf43,r_err_cf43,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='P',label="-43C")
ax2.errorbar(Q_cf45,r_cf45,r_err_cf45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='P',label="-45C")
ax2.errorbar(Q_cf50,r_cf50,r_err_cf50,linestyle='none',color='k',capsize=0,elinewidth=elw,
             marker='P',label="-50C")
ax2.set_xlabel("Qseitz [keV]")
ax2.set_ylabel("Rate [mHz]")
ax2.set_title("Cf-252",loc='right')
ax2.legend(loc='lower right')
ax2.grid()

# Photoneutrons
ax3 = fig.add_subplot(223)
ax3.plot(linex, flaty, linewidth=3, label='constant fit to background')
ax3.plot(linex,flaty1,'k--',label=r'$1\sigma$ region')
ax3.plot(linex,flaty2,'k--')
ax3.fill_between(linex,flaty1,flaty2,facecolor='gray',alpha=0.1)
ax3.errorbar(Q_bibe41,r_bibe41,r_err_bibe41,linestyle='none',color='m',capsize=0,elinewidth=elw,
             marker='s',label = "Bi-Be -41C")
ax3.errorbar(Q_BiBe45,r_BiBe45,r_err_BiBe45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='s',label = "Bi-Be -45C")
ax3.errorbar(Q_bibe43a,r_BiBe43a,r_err_BiBe43a,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='s',fillstyle='none',label = "Bi-Be -43C, 25 psia only")
ax3.errorbar(Q_bibe43b,r_BiBe43b,r_err_BiBe43b,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='s',label = "Bi-Be -43C, scan")

ax3.errorbar(Q_YBeO45,r_YBeO45,r_err_YBeO45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='D',label="Y-Be -45C")
ax3.errorbar(Q_ybe43,r_YBeO43,r_err_YBeO43,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='D',label="Y-Be -43C")
ax3.set_xlabel("Qseitz [keV]")
ax3.set_ylabel("Rate [mHz]")
ax3.set_title("Photoneutron Sources",loc='right')
ax3.legend(loc='lower right')
ax3.grid()

# Gammas
ax4 = fig.add_subplot(224)
ax4.plot(linex, flaty, linewidth=3, label='constant fit to background')
ax4.plot(linex,flaty1,'k--',label=r'$1\sigma$ region for background')
ax4.plot(linex,flaty2,'k--')
ax4.fill_between(linex,flaty1,flaty2,facecolor='gray',alpha=0.1)
ax4.errorbar(Q_bial41,r_bial41,r_err_bial41,linestyle='none',color='m',capsize=0,elinewidth=elw,
             marker='<',markersize=10,label = "Bi-Al -41C")
ax4.errorbar(Q_BiAl45,r_BiAl45,r_err_BiAl45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='<',markersize=10,label = "Bi-Al -45C")


ax4.errorbar(Q_yal45,r_yal45,r_err_yal45,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='^',markersize=10,label = "Y-Al -45C")
ax4.errorbar(Q_yal41,r_yal41,r_err_yal41,linestyle='none',color='m',capsize=0,elinewidth=elw,
             marker='^',markersize=10,label = "Y-Al -41C")
ax4.errorbar(Q_yal38,r_yal38,r_err_yal38,linestyle='none',color='r',capsize=0,elinewidth=elw,
             marker='^',markersize=10,label = "Y-Al -38C")

ax4.errorbar(Q_co45b,r_co45b,r_err_co45b,linestyle='none',color='b',capsize=0,elinewidth=elw,
             marker='v',markersize=10,label = "Co-57 -45C")
ax4.errorbar(Q_co43,r_co43,r_err_co43,linestyle='none',color='g',capsize=0,elinewidth=elw,
             marker='v',markersize=10,label = "Co-57 -43C")
ax4.errorbar(Q_co38,r_co38,r_err_co38,linestyle='none',color='r',capsize=0,elinewidth=elw,
             marker='v',markersize=10,label = "Co-57 -38C")
ax4.set_xlabel("Qseitz [keV]")
ax4.set_ylabel("Rate [mHz]")
ax4.set_title("Gamma Sources",loc='right')
ax4.legend(loc='lower right')
ax4.grid()

plt.show()