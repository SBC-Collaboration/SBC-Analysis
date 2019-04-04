#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:22:37 2019

@author: bressler
"""


import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy
from random import randrange
import matplotlib.dates as dates

def getSBCPMTBaseline(runpath,event,plotbool):
    """
    Returns the list [mu,sigma] for a gaussian fit to the 
    baseline distribution of a SBC event
    """
    ev_baselines=[]
    e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
    tr = e["PMTtraces"]
    trac = tr["traces"]

    for i in range(len(trac)):
        trace = trac[i][0]
        base = np.mean(trace[0:100])
        ev_baselines.append(base)
    def gaus(x,mu,sigma):
            return scipy.stats.norm.pdf(x,mu,sigma)
    
    vals, bins, _ = plt.hist(ev_baselines,int(np.floor(np.sqrt(len(ev_baselines)))),normed=True)
    b = bins[:len(bins)-1]
    params,params_cov = scipy.optimize.curve_fit(gaus,b,vals,p0=[-16,1])
    if plotbool:
        plt.figure()
        x=np.arange(min(ev_baselines),max(ev_baselines),0.1)
        plt.plot(x,[gaus(a,params[0],params[1]) for a in x],c='r',linewidth=3)
        plt.xlabel('trace baselines (mV)',fontsize=15)
        plt.show
    return params
    
def main():
    
    runs902 = ["20170706_4","20170719_0","20170928_11","20171002_3","20171010_10"]
    xAxisDates = [dates.datestr2num((x.split('_'))[0]) for x in runs902]
    baselines=[]
    errors=[]
    for run in runs902:
        run_baselines = []
        run_err = []
        runpath = "/bluearc/storage/SBC-17-data/"+run
        events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
        for e in events:
            m,sigma = getSBCPMTBaseline(runpath,e,False)
            run_baselines.append(m)
            run_err.append(sigma)
        baselines.append(np.mean(run_baselines))
        errors.append(np.mean(run_err))
    plt.figure()
    plt.grid(True)
    plt.plot_date(xAxisDates,baselines,markersize=12)
    plt.xlim([dates.datestr2num('20170701'), dates.datestr2num('20171020')])
    plt.ylabel('PMT Baseline (mV)',fontsize=18)
    #plt.errorbar(np.arange(0,len(runs902)),baselines,errors)
    plt.show
    
if __name__=="__main__":
    main()