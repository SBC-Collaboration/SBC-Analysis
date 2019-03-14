# -*- coding: utf-8 -*-
import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy

runpath = "/bluearc/storage/SBC-17-data/20170706_4/"
events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]

e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,0)

#d=sbc.AnalysisModules.PMTfastDAQalignment.PMTandFastDAQalignment(e)
#print(d.keys())
tr = e["PMTtraces"]
trac = tr["traces"]
t0 = tr["t0"]
dt = tr["dt"]
pmttrig = e["fastDAQ"]["PMTtrig"]
print(tr.keys())

x=np.arange(len(othertrace))
plt.figure()
plt.plot(x,othertrace)
plt.plot(x,trace)
plt.show

V = []
VwithNIM = []
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
plt.hist(np.asarray(V),80,color='r',histtype = 'step')
plt.hist(np.asarray(VwithNIM),80,color='b',histtype='step')
plt.title('RunType 902: 20170706_4, event 0')
plt.xlabel('V max')
plt.show

print(len(indices_of_high))

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
