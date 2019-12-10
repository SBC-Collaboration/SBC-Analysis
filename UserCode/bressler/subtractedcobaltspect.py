#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:59:32 2019

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import matplotlib.pyplot as plt
import scipy
from cobalt57spectrumattempt import spectrum

def main():
    bgrun = '20171003_1'
    bgrunpath = "/bluearc/storage/SBC-17-data/"+bgrun+'/'
    sourcerun = '20171010_6'
    sourcerunpath = '/bluearc/storage/SBC-17-data/'+sourcerun+'/'
    bins = np.arange(-30,130)
    bgspect,_,bgtime = spectrum(bgrunpath,forcebins=bins)
    
    bgspect /= bgtime
    

    sourcespect,sbins,sourcetime = spectrum(sourcerunpath,forcebins=bins)
    print(sourcespect)
    sourcespect /= sourcetime

    ibg = scipy.integrate.trapz(bgspect)
    isource = scipy.integrate.trapz(sourcespect)
    
    print("bg trigger rate: "+str(ibg)+" Hz")
    print("source trigger rate: "+str(isource)+" Hz")


    
    plt.figure()
    plt.plot(bins[:len(sourcespect)],sourcespect,'r',linewidth=4,label="source")
    plt.plot(bins[:len(bgspect)],bgspect,'g',linewidth=4,label='background')
    plt.xlabel('Photoelectrons',fontsize=25)
    plt.ylabel('Rate [Hz/phe]',fontsize=25)
    plt.yscale('log')
    plt.legend(fontsize=18)
    plt.show
    
    subtractedSpectrum = sourcespect - bgspect
    
    def gaussian(x,mu,sigma,amplitude):
        return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )
    
    fitx = [x for x in bins if x> 3 and x <20]
    fitspect = []
    for i in range(len(subtractedSpectrum)):
        if bins[:len(sourcespect)][i] in fitx:
            fitspect.append(subtractedSpectrum[i])

    params, params_cov = scipy.optimize.curve_fit(gaussian,fitx,fitspect, p0=[10,1,10])
    fitx = np.arange(20)
    plt.figure()
    plt.plot(bins[:len(subtractedSpectrum)],subtractedSpectrum,linewidth=5)
    plt.plot(fitx,[gaussian(x,params[0],params[1],params[2]) for x in fitx],'r',linewidth=7)
    plt.title('mean='+str(params[0]))
    plt.xlabel('Photoelectrons',fontsize=18)
    plt.ylabel('Excess in Rate during source run [Hz/phe]',fontsize=18)
    plt.show
    

if __name__ == "__main__":
    main()
