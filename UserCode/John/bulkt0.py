
import os
import matplotlib.pyplot as plt

import numpy as np

import gc

from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE


runid = "20170928_0"

events = [18,19,20,
          21,22,23,
          24,25,26]

acoustic_data = RB("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output/{runid}/AcousticAnalysis_{runid}.bin".format(runid=runid))
fig,ax = plt.subplots(nrows=3, ncols=3)
ax = ax.flatten()
for ev in range(len(events)):
    gc.collect()
    evdata = GE("/bluearc/storage/SBC-17-data/{runid}/".format(runid=runid), events[ev], "fastDAQ")
    time = evdata["fastDAQ"]["time"]
    piezo = evdata["fastDAQ"]["Piezo2"]
    ax[ev].plot(time, piezo, color="b", zorder=1)
    t0 = acoustic_data["bubble_t0"][events[ev]][1]
    if not np.isnan(t0):
        plt.axvline(t0, color="m", linewidth=3, zorder=3)
