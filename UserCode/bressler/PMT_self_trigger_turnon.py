#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 08:31:00 2019

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
from SBCPMTBaseline import getSBCPMTBaseline

def getSBCNIMTurnOn(rp,event,plotbool):
    """
    Returns the list [mu,sigma] for a gaussian CDF fit to the 
    NIM Turn-on/Efficiency distribution of a SBC event
    """

    V = []
    VwithNIM = []
    e = sbc.DataHandling.GetSBCEvent.GetEvent(rp,event)
    tr = e["PMTtraces"]
    trac = tr["traces"]

    for i in range(len(trac)):
        trace = trac[i][0]
        othertrace = trac[i][1]
        V.append(np.fabs(min(trace)))

        if min(othertrace) < -30:
            VwithNIM.append(np.fabs(min(trace)))
    
    vvals, bins, _= plt.hist(np.asarray(V),126,color='r',histtype = 'step')
    vnimvals, _, _ = plt.hist(np.asarray(VwithNIM),bins=bins,color='b',histtype='step')

    for i in range(len(trac)):
        trace = trac[i][0]
        vnimvals = vnimvals[vvals>0]
        vvals = vvals[vvals>0]
    
    perc = np.divide(vnimvals,vvals)
    perc[np.isnan(perc)]=float('+inf')
    perc=perc[perc<float('+inf')]
    
    def functn(x,a,b):
        return scipy.stats.norm.cdf(x,a,b)
    
    params, params_cov = scipy.optimize.curve_fit(functn,bins[:len(perc)],perc,p0=[50,1])
    
    if plotbool:
        plt.figure()
        x=np.arange(min(perc),max(perc),0.1)
        plt.plot(x,[functn(a,params[0],params[1]) for a in x],c='r',linewidth=3)
        plt.xlabel('trace baselines (mV)',fontsize=15)
        plt.show
    return params


def main():
    
    runs902 = ["20170706_4","20170719_0","20170928_11","20171002_3","20171010_10"]
    xAxisDates = [dates.datestr2num((x.split('_'))[0]) for x in runs902]
    TVs=[]
    TVs_rel = []
    errors=[]
    for run in runs902:
        run_tv = []
        run_tv_rel = []
        run_err = []
        runpath = "/bluearc/storage/SBC-17-data/"+run
        events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
        for e in events:
            b,s = getSBCPMTBaseline(runpath,e,False)
            m,sigma = getSBCNIMTurnOn(runpath,e,False)
            run_tv.append(-m)
            run_tv_rel.append(-m - b)
            run_err.append(sigma)
        TVs.append(np.mean(run_tv))
        TVs_rel.append(np.mean(run_tv_rel))
        errors.append(np.mean(run_err))
    print(TVs)
    print(TVs_rel)
    plt.figure()
    plt.grid(True)
    plt.plot_date(xAxisDates,TVs_rel,markersize=12,c='r')
    plt.plot_date(xAxisDates,TVs,markersize=12)
    plt.xlim([dates.datestr2num('20170701'), dates.datestr2num('20171020')])
    plt.ylabel('PMT Self-Trigger Turn On (mV)',fontsize=18)
    plt.legend(['With Baseline Subtraction','Without Baseline Subtraction'])
    #plt.errorbar(np.arange(0,len(runs902)),baselines,errors)
    plt.show
    
if __name__=="__main__":
    main()