#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 09:48:29 2020

@author: bressler
"""

from PICOcode.REFPROP import SeitzModel
import numpy as np
import math
import matplotlib.pyplot as plt

def BTM(p,T,fluid):
    pi = math.pi
    m = SeitzModel.SeitzModel(p,T,fluid)
    Qseitz = m.Q
    rc = m.Rc*1e-9
    sigma = m.Sigma
    dsdT = m.dSigmadT
    rhoL = m.Rho_l * 1000
    rhoB = m.Rho_b * 1000
    A = 17
    
    A_Auger = 3 # per k-shell photoabsorption
    B_C3F8I = (1/230) # MeV cm^2/g
    
    if fluid == 'r218' or fluid == 'C3F8' or fluid == 'c3f8':
        B = 0.02703
    elif fluid == 'xenon':
        B = 0.09615
    elif fluid == 'argon':
        B = 0.04166
    elif fluid == 'ethylene':
        B = 0.019
    p = p * 6894.76
    T = T + 273.15
    
    Eion = 4*pi*rc*rc*(sigma - T*dsdT) + (4*pi/3)*rc*rc*rc*p
    
    Eion = Eion * 6.242e12 # MeV
    
    rl = rc*((rhoB/rhoL)**(1/3))*100 # cm
    
    f_ion = (Eion/(rl*rhoL/1000)) # MeV cm^2 g^-1 
    P_ion = (A*np.exp(-B*f_ion))*1000 # events per keV
    
    f_seitz = Qseitz/(rl*rhoL) # MeV cm^2 g^-1
    P_contaminated = (A_Auger*np.exp(-B_C3F8I * f_seitz)) # events per k-shell photoabsorption
    
    return [Qseitz, Eion, f_ion, P_ion, f_seitz, P_contaminated]
    
def main():
    print(BTM(20, 14.5, 'r218'))
    """
    with open('/nashome/b/bressler/sbcoutput/PICOC3F8_25psia_gamma.txt', 'w') as fout:
        C3F8_25_p = np.zeros([50,1])
        for i in range(len(np.linspace(10,25,50))):
            T = np.linspace(10,25,50)[i]
            C3F8_25_p[i] = BTM(25,T,'r218')[3]
            fout.write("%f %e\n"%(BTM(25,T,'r218')[0],BTM(25,T,'r218')[3]))
    print(C3F8_25_p)
    
    with open('/nashome/b/bressler/sbcoutput/DBC_xeSpike_Q_and_f.txt', 'w') as fout:
        T = 19
        for p in np.arange(15,55,0.1):
            btm = BTM(p, T, 'r218')
            fout.write("%f %e %e %e %e %e\n"%(btm[0], btm[1], btm[2], btm[3], btm[4], btm[5]))
    
    
    plt.figure()
    Tarr = np.array([-45,-40,-35,-30,-25,-20])
    #Tarr = np.array([-25,-20,-18])
    Tarrar = np.array([-138, -143,-145,-150,-155])
    Parr = np.arange(15,100)
    Qxe = np.array([np.zeros([1,len(Parr)]) for i in range(len(Tarr))])
    pxe = np.array([np.zeros([1,len(Parr)]) for i in range(len(Tarr))])
    Qar = np.array([np.zeros([1,len(Parr)]) for i in range(len(Tarr))])
    par = np.array([np.zeros([1,len(Parr)]) for i in range(len(Tarr))])
    for i in range(len(Tarr)):
        for j in range(len(Parr)):
            Qxe[i][0][j]=BTM(float(Parr[j]),float(Tarr[i]),'xenon')[0]
            pxe[i][0][j]=BTM(float(Parr[j]),float(Tarr[i]),'xenon')[3]
        plt.plot(Qxe[i][0][:], pxe[i][0][:],'--',label='%f C, xenon'%Tarr[i])
    for i in range(len(Tarrar)):
        for j in range(len(Parr)):
            Qar[i][0][j]=BTM(float(Parr[j]),float(Tarrar[i]),'argon')[0]
            par[i][0][j]=BTM(float(Parr[j]),float(Tarrar[i]),'argon')[3]
        plt.plot(Qar[i][0][:],par[i][0][:],label='%f C, argon'%Tarrar[i])
    plt.yscale('log')
    plt.xscale('log')
    plt.grid('both')
    plt.xlabel('Qseitz [keV]')
    plt.ylabel('probability of bubble nucleation per keV deposited')
    plt.legend()
    plt.show()
    
    print(BTM(145, -138, 'argon'))
    """
if __name__=="__main__":
    main()