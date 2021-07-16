#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:09:35 2019

@author: bressler
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

rho_D = 0.4 #GeV c^-2 cm^-3
v0 = 230 #km/s
v_esc = 600 #km/s
N0 = 6.02e26 #mol^-1 , avagadro's number
A_Ar = 39.948 #u, atomic mass of argon

r = lambda MD,MT: (4*MD*MT) / ((MD+MT)**2)

f = lambda v,vE: np.exp((-(v+vE)**2) / v0**2)



k0 = (np.pi * v0**2)**(3/2)
k = k0 * (sp.special.erf(v_esc/v0) - (2/(np.sqrt(np.pi))) * (v_esc/v0) * np.exp(-(v_esc**2)/v0**2))

def simple_dRdER(A,MD,ER,sigma):
    MT = 0.932*A
    R0 = (540/(A*MD)) * sigma
    E0 = 0.5 * MD * v0**2
    return (R0/E0) * np.exp(-ER/(E0*r(MD,MT)))

Estep = 1
Earr = np.arange(0.1,1e5,Estep)
dR = [simple_dRdER(A_Ar,1,Earr[i],1) * Estep for i in range(len(Earr))]
R = np.cumsum(dR)

plt.figure()
plt.plot(Earr,R)
plt.show
plt.figure()
plt.plot(Earr,dR)
plt.show