# -*- coding: utf-8 -*-
import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy

runpath = "/bluearc/storage/SBC-17-data/20170719_0/"
events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]

V = []
VwithNIM = []
for event in events:
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)

    #d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
    #print(d.keys())
    tr = e["PMTtraces"]
    trac = tr["traces"]
    t0 = tr["t0"]
    dt = tr["dt"]
    pmttrig = e["fastDAQ"]["PMTtrig"]


    indices_of_high = []
    for i in range(len(trac)):
        trace = trac[i][0]
        othertrace = trac[i][1]
        V.append(np.fabs(min(trace)))
    
        if np.fabs(min(trace)) == 128:
            indices_of_high.append(i)
        if min(othertrace) < -30:
            VwithNIM.append(np.fabs(min(trace)))
plt.figure()
vvals, bins, _= plt.hist(np.asarray(V),110,color='r',histtype = 'step')
vnimvals, _, _ = plt.hist(np.asarray(VwithNIM),bins=bins,color='b',histtype='step')
plt.title('RunType 902: 20170719_0')
plt.xlabel('V max')
plt.show

vnimvals = vnimvals[vvals>0]
vvals = vvals[vvals>0]

diff = vvals-vnimvals
perc = np.divide(vnimvals,vvals)
perc[np.isnan(perc)]=float('+inf')
perc=perc[perc<float('+inf')]

def functn(x,a,b):
    return scipy.stats.norm.cdf(x,a,b)

params, params_cov = scipy.optimize.curve_fit(functn,bins[:len(perc)],perc,p0=[50,1])

plt.figure()
plt.scatter(bins[:len(perc)],perc)
plt.plot(bins[:len(perc)],functn(bins[:len(perc)],params[0],params[1]),color='r')
plt.text(60,.75,"mu = "+str(params[0]),fontsize=15)
plt.text(60,.5,"sigma = "+str(params[1]),fontsize=15)
plt.xlabel('V max')
plt.ylabel('efficiency')
plt.show()
"""
for j in range(len(indices_of_high)):
    trace = np.fabs(trac[indices_of_high[j]][0])
    baseline = np.mean(trace[0:100])
    trace = trace - baseline
    pk_ind = scipy.signal.find_peaks(trace,5)
    pk_vals = [trace[k] for k in pk_ind[0]]
    plt.figure()
    plt.hold(True)
    x=np.arange(len(trace))
    plt.plot(x,trace)
    plt.scatter(pk_ind[0],pk_vals,s=50,c="r")
    plt.show
"""