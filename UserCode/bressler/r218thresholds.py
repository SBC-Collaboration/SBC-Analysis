#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 10:29:20 2019

@author: bressler
"""

from PICOcode.REFPROP.SeitzModel import SeitzModel
import numpy as np
import matplotlib.pyplot as plt

T = np.arange(-153,-138,1)
#T=[-143]
print(T)
P= np.arange(15,382,2)

QperT = [[] for i in range(len(T))]
RperT = [[] for i in range(len(T))]

for i in range(len(T)):
    print(T[i])
    for j in range(len(P)):
        m=SeitzModel(float(P[j]),float(T[i]),'argon')
        QperT[i].append(m.Q[0])
        RperT[i].append(m.Rc[0])
        
plt.figure()
for i in range(len(T)):
    plt.plot(QperT[i],RperT[i],label="T=%f K"%(T[i]+273.15))
plt.xlim([0.01,0.5])
plt.ylim([0,15])
plt.xlabel("Seitz Threshold [keV]")
plt.ylabel("Critical Radius [nm]")
plt.legend()
plt.show()