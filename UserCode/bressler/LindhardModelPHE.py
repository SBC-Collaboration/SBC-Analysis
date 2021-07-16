#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 14:15:41 2020

@author: bressler
"""

import numpy as np
import matplotlib.pyplot as plt

efficiency = 0.0005

def calculateXenonNeNph(E):
    """ Give E in keV"""
    k = 0.1735
    W = 13.7
    eps = 11.5 * (54**(-7/3)) * E
    g = 3 * (np.sign(eps)*np.abs(eps)**0.15) + 0.7 * (np.sign(eps)*np.abs(eps)**0.6) + eps
    L = k * g / (1 + k * g)
    
    NeNph = (L*E*1000)/W
    return NeNph

def pheXenonChamber(E):
    return efficiency * calculateXenonNeNph(E)

def EnrXenonChamber(phe):
    energyarray = np.arange(0,1e6,0.1)
    return np.interp(phe, pheXenonChamber(energyarray), energyarray)

def main():
    Enr = np.arange(0.1, 1e6, 0.1)
    NeNph = np.zeros([len(Enr),1])
    print(EnrXenonChamber(1))
    print(EnrXenonChamber(512.5))
    print(pheXenonChamber(15000))
    
    for i in range(len(Enr)):
        E = Enr[i]
        NeNph[i] = calculateXenonNeNph(E)
    
    plt.figure()
    plt.plot(Enr, [efficiency*N for N in NeNph])
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.xlabel(r'$E_{nr}$ [keV]', fontsize=18)
    plt.ylabel('Expected Photoelectrons', fontsize=18)
    plt.title(r'$\epsilon=%f$'%efficiency, fontsize=18)
    plt.show()

if __name__=="__main__":
    main()