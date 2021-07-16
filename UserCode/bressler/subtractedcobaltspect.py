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
    bgrun = '20171002_8'
    datapath =  "/bluearc/storage/SBC-17-data/"
    bgrunpath = datapath + bgrun+'/'
    sourcerun = '20171002_4'
    sourcerunpath = datapath + sourcerun+'/'
    bins = np.arange(-1,300,0.5)
    bgspect,_,bgtime = spectrum(datapath, bgrun,forcebins=bins)
    
    bgspect /= bgtime
    binwidths = [bins[i+1]-bins[i] for i in range(len(bins)-1)]
    bincenters = [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1)]
    bgdiffrate = [bgspect[i]/binwidths[i] for i in range(len(bgspect))]
    print(bgdiffrate)
    print(bincenters)
    sourcespect,sbins,sourcetime = spectrum(datapath, sourcerun,forcebins=bins)
    sbinwidths = [sbins[i+1]-sbins[i] for i in range(len(sbins)-1)]
    sourcespect /= sourcetime
    sourcediffrate = [sourcespect[i]/sbinwidths[i] for i in range(len(sourcespect))]
    sbincenters = [0.5*(sbins[i+1]+sbins[i]) for i in range(len(sbins)-1)]

    ibg = scipy.integrate.trapz(bgspect)
    isource = scipy.integrate.trapz(sourcespect)
    
    print("bg trigger rate: "+str(ibg)+" Hz")
    print("source trigger rate: "+str(isource)+" Hz")


    print(len(bgspect))
    print(len(sourcespect))
    
    plt.figure()
    plt.plot(sbincenters,sourcespect,'r',linewidth=4,label=sourcerun)
    plt.plot(bincenters,bgspect,'g',linewidth=4,label=bgrun)
    plt.xlabel('Photoelectrons',fontsize=25)
    plt.ylabel('Rate [Hz/bin]',fontsize=25)
    plt.yscale('log')
    plt.legend(fontsize=18)
    plt.show
    
    plt.figure()
    plt.plot(sbincenters,sourcediffrate,'r',linewidth=4,label=sourcerun)
    plt.plot(bincenters,bgdiffrate,'g',linewidth=4,label=bgrun)
    plt.xlabel('Photoelectrons',fontsize=25)
    plt.ylabel('Diff. Rate [Hz/phe]',fontsize=25)
    plt.yscale('log')
    plt.legend(fontsize=18)
    plt.show    
    
    subtractedSpectrum = [sourcediffrate[i] - bgdiffrate[i] for i in range(len(sourcediffrate))]
    
    
    def gaussian(x,mu,sigma,amplitude):
        return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )
    
    def addedGaussians(x, mu1, mu2, sigma1, sigma2, amplitude1, amplitude2):
        return amplitude1 * np.exp(-((x - mu1) /(np.sqrt(2)* sigma1))**2 ) + amplitude2 * np.exp(-((x - mu2) /(np.sqrt(2)* sigma2))**2 )
    
    highlighthumpx = [x for x in bincenters if x >40 and x< 100]
    highlightspect = []
    
    fitx = [x for x in bincenters if x> 0 and x <20]
    plotx = np.arange(0,max(fitx),0.01)
    fitspect = []
    for i in range(len(subtractedSpectrum)):
        if sbincenters[i] in fitx:
            fitspect.append(subtractedSpectrum[i])
        if sbincenters[i] in highlighthumpx:
            highlightspect.append(subtractedSpectrum[i])
            
    #highparams, highparams_cov = scipy.optimize.curve_fit(gaussian, highlighthumpx, highlightspect, p0=[110,70,0.5])

    params, params_cov = scipy.optimize.curve_fit(addedGaussians,fitx,fitspect, p0=[1,2,1,10,30,10])
    
    plt.figure()
    plt.plot(sbincenters,subtractedSpectrum,linewidth=1,label='data')
    plt.plot(plotx,[addedGaussians(x,params[0],params[1],params[2],params[3],params[4],params[5]) for x in plotx],'r:',
                   linewidth=3,label=r'two gaussian fit, $\mu_1=%.2f$, $\mu_2 = %.2f$'%(params[0], params[1]))
    #plt.plot(highlighthumpx, [gaussian(x, highparams[0], highparams[1], highparams[2]) for x in highlighthumpx],'--b',linewidth=4,label=r'$\mu$=%f'%highparams[0])
    plt.title(sourcerun+'-'+bgrun+r', 122 keV $\mu=%.2f$'%(params[1])+r", $\sigma=%.2f$ "%(params[3]))
    plt.xlabel('Photoelectrons',fontsize=18)
    plt.ylabel('Excess in Rate during source run [Hz/phe]',fontsize=18)
    plt.yscale('log')
    plt.legend()
    plt.show
    print("122 keV = %f phe"%params[1])

if __name__ == "__main__":
    main()
