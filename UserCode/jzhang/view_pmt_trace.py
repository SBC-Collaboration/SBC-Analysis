import numpy as np
import os
import re
import matplotlib.pyplot as plt

import SBCcode as sbc


def loaddata(runid, ev):
    datadir = '/mnt/XENON_DAQ/SBC-17-data'
    loadlist = ['fastDAQ', 'PMTtraces']
    # thisevent = sbc.get_event(os.path.join(datadir, runid), ev, *loadlist)
    # loadlist = ['images','PMTtraces']
    thisevent = sbc.get_event(os.path.join(datadir, runid), ev, *loadlist)
    return thisevent


def plottrace(thisevent, channel=0, ni=1, ne=10):
    plt.figure()
    for i in range(ni, ne + 1):
        plt.clf()
        plt.plot(thisevent['PMTtraces']['traces'][i, channel, :])
        plt.draw()
        plt.waitforbuttonpress()

    # pdb.set_trace()


def plot_fastDAQ(ev):
    if ev['fastDAQ']['loaded']:
        # ['Dytran', 'multiboards', 'Piezo1', 'caldata', 'CAMgate', 'PMTtrig', 'bindata', 'time', 'loaded', 'VetoCoinc']
        keys = ['Dytran', 'Piezo1', 'Piezo2', 'CAMgate', 'PMTtrig', 'time', 'VetoCoinc']
        keys = set(keys).intersection(ev['fastDAQ'].keys())
        # plt.ion()
        for key in keys:
            print(key)
            plt.figure()
            plt.plot(ev['fastDAQ'][key], label=key)
            plt.legend()
            plt.show()

# loaddata('20170630_7',0)
