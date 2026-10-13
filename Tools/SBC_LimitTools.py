#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:44:30 2026

@author: cdahl
"""

import numpy as np
import scipy.optimize as opt

def PoissonChisq(n, nu):
    ''' Calculates Poisson chi-square '''
    # see Baker and Cousins for reference
    chisq_0 = 2 * (nu - n)
    chisq_1 = -2 * n * np.log(nu / n)
    chisq_1[np.isnan(chisq_1)] = 0
    return chisq_0 + chisq_1

def PoissonChi(n, nu):
    ''' Calculates signed sqrt of Poisson chi-square '''
    return np.sqrt(PoissonChisq(n, nu)) * np.sign(nu - n)

def PoissonChisqWBkg(n, nu, k, r):
    ''' Calculates Poisson chi-square with bkg dataset
        n = # of events observed in "signal" dataset
        nu = # of signal events expected
        k = # of events observed in "background" dataset
        r = exposure of "signal" dataset relative to "background" '''
        
    b = n + k - nu * ((1 + r) / r)
    a = 1 + r
    mu = (b + np.sqrt(b * b + 4 * a * k * (nu / r))) / (2 * a)
    mu[r == 0] = k[r == 0]  # special case -- known background
    chisq_0 = 2 * (nu + mu + mu * r - n - k)
    chisq_1 = -2 * n * np.log((nu + mu * r) / n)
    chisq_1[np.isnan(chisq_1)] = 0
    chisq_2 = -2 * k * np.log(mu / k)
    chisq_2[np.isnan(chisq_2)] = 0
    return chisq_0 + chisq_1 + chisq_2

def PoissonChiWBkg(n, nu, k, r):
    ''' Calculates signed sqrt of Poisson chi-square w bkg '''
    return np.sqrt(PoissonChisqWBkg(n, nu, k, r)) * np.sign(nu - n + k * r)

def PoissonChisqWBkg_UL(n, nu, k, r):
    ''' Calculates Poisson chi-square with bkg dataset
        n = # of events observed in "signal" dataset
        nu = # of signal events expected
        k = # of events observed in "background" dataset
        r = exposure of "signal" dataset relative to "background" 
        BUT -- returns zero if nu < (n - k*r) '''
        
    b = n + k - nu * ((1 + r) / r)
    a = 1 + r
    mu = (b + np.sqrt(b * b + 4 * a * k * (nu / r))) / (2 * a)
    mu[r == 0] = k[r == 0]  # special case -- known background
    chisq_0 = 2 * (nu + mu + mu * r - n - k)
    chisq_1 = -2 * n * np.log((nu + mu * r) / n)
    chisq_1[np.isnan(chisq_1)] = 0
    chisq_2 = -2 * k * np.log(mu / k)
    chisq_2[np.isnan(chisq_2)] = 0
    chisq_tot = chisq_0 + chisq_1 + chisq_2
    return chisq_tot * (nu > (n - k*r))

def FindUL(n, nu_ref, k, r, target_chisq):
    ''' Calculates upper limit given
        n = # of events observed in "signal" dataset
        nu_ref = # of signal events expected for some reference cross-section
        k = # of events observed in "background" dataset
        r = exposure of "signal" dataset relative to "background" 
             --- all of n, nu_ref, k, r, should have shape(N,) ---
        target_chisq = the chisq (as calcualted by PoissonChisqWBkg_UL)
               associated with this upper limit (must be calcualted separately)
        with OUTPUT:
        sol = upper limit on scale factor for nu_ref '''
        
    initial_guess = np.max((n + .5 + np.sqrt(.25+n)) / nu_ref)
    # initial_guess is at least +1 sigma for every point before bkg subtraction
    sol = opt.fsolve(lambda x: np.sum(PoissonChisqWBkg_UL(n, x*nu_ref, k, r) - target_chisq),
                     initial_guess)
    return sol[0]
    
def GenerateChisqDist(nu_sig, nu_bkg, r, ntoys):
    ''' Use PoissonChisqWBkg_UL to calculate a bunch of chisquares
        given a signal and background model (assuming no untagged bkg) --- 
        This is meant to help find the target_chisq for the FindUL function '''
    toy_ones = np.ones((1,ntoys), dtype=np.float64)    
    n = np.random.poisson((nu_sig + nu_bkg*r)[:,None]*toy_ones)
    k = np.random.poisson(nu_bkg[:,None]*toy_ones)
    chisq = np.zeros(n.shape, dtype=np.float64)
    for i_t in np.arange(ntoys, dtype=np.intp):
        chisq[:,i_t] = PoissonChisqWBkg_UL(n[:,i_t], nu_sig, k[:,i_t], r)
    return np.sum(chisq,axis=0)

