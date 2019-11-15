#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 10:24:55 2019

@author: bressler
"""

import SBCcode as sbc
from os import listdir
from os.path import isfile,join
import numpy as np
import matplotlib.pyplot as plt
import scipy


def NIM_efficiency_and_plot(V,VwithNIM,title_string):
    plt.figure()
    vvals, bins, _= plt.hist(np.asarray(V),126,color='r',histtype = 'step')
    vnimvals, _, _ = plt.hist(np.asarray(VwithNIM),bins=bins,color='b',histtype='step')

    plt.xlabel('V max')
    plt.title(title_string)
    plt.show
    
    vnimvals = vnimvals[vvals>0]
    vvals = vvals[vvals>0]
    
    perc = np.divide(vnimvals,vvals)
    perc[np.isnan(perc)]=float('+inf')
    perc=perc[perc<float('+inf')]
    
    def functn(x,a,b):
        return scipy.stats.norm.cdf(x,a,b)
    
    params, params_cov = scipy.optimize.curve_fit(functn,bins[:len(perc)],perc,p0=[50,1])
    
    plt.figure()
    plt.scatter(bins[:len(perc)],perc)
    plt.plot(bins[:len(perc)],functn(bins[:len(perc)],params[0],params[1]),color='r')
    plt.text(40,.75,"mu = "+str(params[0]),fontsize=15)
    plt.text(40,.5,"sigma = "+str(params[1]),fontsize=15)
    plt.xlabel('V max')
    plt.ylabel('efficiency')
    plt.title(title_string)
    plt.show()
    
    
def main():
    runpath = "/bluearc/storage/SBC-17-data/20170719_0/"
    events = [evnt for evnt in listdir(runpath) if not isfile(join(runpath,evnt))]
    
    V = []
    VwithNIM = []
    for event in events:
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runpath,event)
        tr = e["PMTtraces"]
        trac = tr["traces"]
    
        for i in range(len(trac)):
            trace = trac[i][0]
            othertrace = trac[i][1]
            V.append(np.fabs(min(trace)))
    
            if min(othertrace) < -30:
                VwithNIM.append(np.fabs(min(trace)))
    NIM_efficiency_and_plot(V,VwithNIM,'title')
    
if __name__=="__main__":
    main()
    

