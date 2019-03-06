# -*- coding: utf-8 -*-
import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy

runpath = "/bluearc/storage/SBC-17-data/20170718_4/"
events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]

e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,0)
tr = e["PMTtraces"]
trac = tr["traces"]
t0 = tr["t0"]
dt = tr["dt"]
print(tr.keys())
print(t0[0])

V = []
indices_of_high = []
for i in range(len(trac)):
    trace = trac[i][0]
    V.append(np.fabs(min(trace)))
    
    if np.fabs(min(trace)) == 128:
        indices_of_high.append(i)
plt.figure()
plt.hist(np.asarray(V),100)
plt.show

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
