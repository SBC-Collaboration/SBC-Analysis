#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 09:30:30 2020

@author: bressler
"""
from PICOcode.REFPROP.SeitzModel import SeitzModel
import numpy as np
import matplotlib.pyplot as plt
from baxterelectronrecoilmodel import BTM

with open('/coupp/data/home/coupp/users/bressler/output/xenonspikeoutput.txt','r') as xenonfile:
    d = xenonfile.readlines()

data = {}
for line in d:
    elements = line.split()
    #print([float(e) for e in elements])
    
    if int(elements[1]) == 0:
        if float(elements[5]) <=18:
            data["Background 17C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
        else:
             data["Background 19C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
           
    elif int(elements[1]) == 137:
        if float(elements[5]) <=18:
            data["Cs-137 17C position %d %.1f ppm %.1f psia"%(int(elements[2]), float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
        else:
            data["Cs-137 19C position %d %.1f ppm %.1f psia"%(int(elements[2]), float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]

    elif int(elements[1]) == 152:
        if float(elements[5]) <=18:
            data["Eu-152 17C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
        else:
            data["Eu-152 19C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
    
    elif int(elements[1]) == 244:
        if float(elements[5]) <=18:
            data["Cm-244 17C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
        else:
            data["Cm-244 19C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
            
    elif int(elements[1]) == 271:
        if float(elements[5]) <=18:
            data["Eu-152 + tin 17C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
        else:
            data["Eu-152 + tin 19C %.1f ppm %.1f psia"%(float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
            
    elif int(elements[1]) == 256:
        if float(elements[5]) <=18:
            data["Cs-137 + tin 17C position %d %.1f ppm %.1f psia"%(int(elements[2]), float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]
        else:
            data["Cs-137 + tin 19C position %d %.1f ppm %.1f psia"%(int(elements[2]), float(elements[0]), float(elements[3]))] = [float(elements[i]) for i in range(len(elements))]

subtractedData = {}
subtractedxenonfile = open('/coupp/data/home/coupp/users/bressler/output/subtractedxenonoutput.txt','w')
for k in data.keys():
    datum = data[k]
    #print(datum)
    if datum[1] != 0.0:
        if datum[5] <= 18:
            if datum[0] == 31.7 and datum[3] == 39.3:
                continue
            else:
                associatedbg = data["Background 17C %.1f ppm %.1f psia"%(float(datum[0]), float(datum[3]))]
        else:
            associatedbg = data["Background 19C %.1f ppm %.1f psia"%(float(datum[0]), float(datum[3]))]
        bgsubrate = datum[7] - associatedbg[7]
        bgsubrateerror = np.sqrt(datum[8]**2 + associatedbg[8]**2)
        #print(bgsubrate)
        [Qseitz, Eion, f_ion, P_ion, f_seitz, P_contaminated]= BTM(datum[4]-1.3, datum[5], 'r218')
        subtractedData[k] = [datum[0], datum[1], datum[2], datum[3], datum[4], datum[5], Qseitz, f_seitz, bgsubrate, bgsubrateerror]
        subtractedxenonfile.write(
                '%.1f %d %d %.1f %.1f %.1f %.1f %.2f %.2f %.2f\n'%(datum[0], 
                    datum[1], datum[2], datum[3], datum[4], datum[5], Qseitz, f_seitz, bgsubrate, bgsubrateerror))
subtractedxenonfile.close()