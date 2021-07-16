#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 14:19:36 2021

@author: bressler
"""

from PICOcode.REFPROP.SeitzModel import SeitzModel
import numpy as np
import matplotlib.pyplot as plt
from baxterelectronrecoilmodel import BTM

with open('/coupp/data/home/coupp/users/bressler/output/argonspikeoutput.txt','r') as argonfile:
    d = argonfile.readlines()

data = {}
for line in d:
    elements = line.split()
    #print([float(e) for e in elements])
    
    if int(elements[1]) == 0:
        data["Background 20C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
           
    elif int(elements[1]) == 137:
        data["Cs-137 20 position %d %.1f ppm %.1f psia"%(int(elements[2]), float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
  
    elif int(elements[1]) == 244:

        data["Cm-244 20 %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]

            

subtractedData = {}
subtractedargonfile = open('/coupp/data/home/coupp/users/bressler/output/subtractedargonoutput.txt','w')
for k in data.keys():
    datum = data[k]
    #print(datum)
    if datum[1] != 0.0:
        associatedbg = data["Background 20C %.1f ppm %.1f psia"%(float(datum[0]), float(datum[3]))]
        bgsubrate = datum[7] - associatedbg[7]
        bgsubrateerror = np.sqrt(datum[8]**2 + associatedbg[8]**2)
        #print(bgsubrate)
        [Qseitz, Eion, f_ion, P_ion, f_seitz, P_contaminated]= BTM(datum[4]-1.3, datum[5], 'r218')
        print("P=%.2f"%datum[4])
        print("T=%.2f"%datum[5])
        print("Qseitz=%.3f"%Qseitz)
        subtractedData[k] = [datum[0], datum[1], datum[2], datum[3], datum[4], datum[5], Qseitz, f_seitz, bgsubrate, bgsubrateerror]
        subtractedargonfile.write(
                '%.1f %d %d %.1f %.1f %.1f %.2f %.2f %.2f %.2f\n'%(datum[0], 
                    datum[1], datum[2], datum[3], datum[4]-1.3, datum[5], Qseitz, f_seitz, bgsubrate, bgsubrateerror))
subtractedargonfile.close()