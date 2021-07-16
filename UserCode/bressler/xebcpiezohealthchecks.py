#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:36:12 2021

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import os
import runlistscatalogue as rlc
import matplotlib.pyplot as plt
import gc


runs = rlc.cfJuly6to11
datadir = '/bluearc/storage/SBC-17-data/'

rise_times_1 = np.zeros(101*len(runs))
rise_times_2 = np.zeros(101*len(runs))
dpk = np.zeros(101*len(runs))

dt = np.zeros(101*len(runs))
rms1 = np.zeros(101*len(runs))
rms2 = np.zeros(101*len(runs))
pretriggerE1 = np.zeros(101*len(runs))
triggerE1 = np.zeros(101*len(runs))
pretriggerE2 = np.zeros(101*len(runs))
triggerE2 = np.zeros(101*len(runs))

pretriggertwinindex = 0
triggerwinindex = 2

frequencybandtouse = 2

totalevents = 0
for run in runs:
    print(run)
    runrawpath = datadir+run
    runreconpath = '/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/'%run
    
    acousticfilename = runreconpath+'AcousticAnalysis_%s.bin'%run
    acousticanalysisdata = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
    #print(acousticanalysisdata.keys())
    

    events = [evnt for evnt in os.listdir(runrawpath) if not os.path.isfile(os.path.join(runrawpath, evnt))]
    Nevents = len(events)
    
    for i in range(Nevents):
        if i+1 > len(acousticanalysisdata['bubble_t0'][:,0]):
            print('incomplete event')
            continue
        e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath, i, 'fastDAQ')
        p1 = e['fastDAQ']['Piezo1']
        p1_noisewindow = p1[:100]
        p1noisesquare = [x**2 for x in p1_noisewindow]
        p1_pk = acousticanalysisdata['peak_t0'][i, 0]
        
        p2 = e['fastDAQ']['Piezo2']
        p2_noisewindow = p2[:100]
        p2noisesquare = [x**2 for x in p2_noisewindow]
        p2_pk = acousticanalysisdata['peak_t0'][i, 1]
        
        at0_1 = acousticanalysisdata['bubble_t0'][i,0]
        at0_2 = acousticanalysisdata['bubble_t0'][i,1]
        
        pte1 = acousticanalysisdata['piezoE'][i, 0, pretriggertwinindex, frequencybandtouse]
        pte2 = acousticanalysisdata['piezoE'][i, 1, pretriggertwinindex, frequencybandtouse]
        
        twe1 = acousticanalysisdata['piezoE'][i, 0, triggerwinindex, frequencybandtouse]
        twe2 = acousticanalysisdata['piezoE'][i, 1, triggerwinindex, frequencybandtouse]
        
        #if twe2 < 1e6:
        #    print('%d has trig window energy p2 less than 1e6'%i)

        delta_at0 = abs(at0_2-at0_1)
        delta_pk = abs(p1_pk - p2_pk)
        rms_1 = np.sqrt(np.mean(p1noisesquare))
        rms_2 = np.sqrt(np.mean(p2noisesquare))
        
        dt[totalevents] = delta_at0
        rms1[totalevents] = rms_1
        rms2[totalevents] = rms_2
        pretriggerE1[totalevents] = pte1
        pretriggerE2[totalevents] = pte2
        
        triggerE1[totalevents] = twe1
        triggerE2[totalevents] = twe2
        
        dpk[totalevents] = delta_pk
        rise_times_1[totalevents] = p1_pk - at0_1
        rise_times_2[totalevents] = p2_pk - at0_2

        totalevents += 1
        gc.collect()

dt = dt[:totalevents]
rms1 = rms1[:totalevents]    
rms2 = rms2[:totalevents]
pretriggerE1 = pretriggerE1[:totalevents]
pretriggerE2 = pretriggerE2[:totalevents]

triggerE1 = triggerE1[:totalevents]
triggerE2 = triggerE2[:totalevents]

dpk = dpk[:totalevents]
rise_times_1 = rise_times_1[:totalevents]
rise_times_2 = rise_times_2[:totalevents]


plt.figure()
plt.scatter(rms1, dt, label='piezo 1')
plt.scatter(rms2, dt, label='piezo 2')
plt.ylabel(r'$|\Delta t_0|$ [s]')
plt.xlabel('rms of piezo noise')
plt.legend()
plt.show()

plt.figure()
plt.scatter(rise_times_2, dt)
plt.xlabel('piezo 2 rise time [s]', fontsize=18)
plt.ylabel(r'$|\Delta t_0|$ [s]', fontsize=18)
plt.grid()
plt.show()

plt.figure()
plt.scatter(rise_times_1, rise_times_2)
plt.xlabel('rise time, piezo 1 [s]')
plt.ylabel('rise time, piezo 2 [s]')
plt.show()

plt.figure()
plt.scatter(dpk, dt)
plt.xlabel(r'$|\Delta t_{rise}|$ [s]', fontsize=18)
plt.ylabel(r'$|\Delta t_0|$ [s]', fontsize=18)
plt.grid()
plt.show()


plt.figure()
plt.scatter(rms1, rms2)
plt.xlabel('rms 1')
plt.ylabel('rms 2')
plt.show()

plt.figure()
plt.scatter(pretriggerE1, dt, label = 'piezo 1')
plt.scatter(pretriggerE2, dt, label = 'piezo 2')
plt.legend(fontsize=15)
plt.xlabel('pretrigger energy', fontsize=18)
plt.ylabel(r'$|\Delta t_0|$ [s]', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-6, 1])
plt.show()

plt.figure()
plt.scatter(triggerE1, dt, label = 'piezo 1')
plt.scatter(triggerE2, dt, label = 'piezo 2')
plt.legend(fontsize=15)
plt.xlabel('trigger window energy', fontsize=18)
plt.ylabel(r'$|\Delta t_0|$ [s]', fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim([1e-6, 1])
plt.show()


with open("/nashome/b/bressler/sbcoutput/cfJuly6to11_merged.txt","r") as fin:
    data = fin.readlines()


headers = data[0].split()
print(headers)
xind = headers.index("x")
yind = headers.index("y")
zind = headers.index("z")
at00ind = headers.index("at0_0")
at01ind = headers.index("at0_1")
pmtt0ind = headers.index("PMTt0")
lagind = headers.index("lag")
spectind = headers.index("PMTphe")
blockedind = headers.index('isBlocked')
x=[]
y=[]
z=[]
lag = []
delta_at0 = []
triggerwindowenergyp2 = []
Vbaseline1 = []
Vbaseline2 = []

dpk = []
rise_times_1 = []
rise_times_2 = []

i=0
for line in data:
    if i>0:
        split_line = line.split()
        if float(split_line[blockedind])>0.9:
            print(split_line[blockedind])
        if (not np.isnan(float(split_line[xind]))):
            
            run = split_line[0]
            evn = int(split_line[1])
            e = sbc.DataHandling.GetSBCEvent.GetEvent('/bluearc/storage/SBC-17-data/%s'%run, evn, 'fastDAQ', 'slowDAQ')
            indices_back = 30
            trig_pressure = e["slowDAQ"]['PT6'][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
            if trig_pressure > 25.5:
                continue
            p1 = e['fastDAQ']['Piezo1']
            
            p2 = e['fastDAQ']['Piezo2']
            x.append(float(split_line[xind]))
            y.append(float(split_line[yind]))
            z.append(float(split_line[zind]))
            lag.append(float(split_line[lagind]))
            delta_at0.append(abs((float(split_line[at00ind])-float(split_line[at01ind]))))
            Vbaseline1.append(np.mean(p1[:1000]))
            Vbaseline2.append(np.mean(p2[:1000]))
            runreconpath = '/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/'%run
            acousticfilename = runreconpath+'AcousticAnalysis_%s.bin'%run
            acousticanalysisdata = sbc.DataHandling.ReadBinary.ReadBlock(acousticfilename)
            twe2 = acousticanalysisdata['piezoE'][evn, 1, triggerwinindex, frequencybandtouse]
            if twe2 < 1e6:
                print('run %s event %d has bubble and trig window energy p2 less than 1e6'%(run, evn))
            triggerwindowenergyp2.append(twe2)
            p1_pk = acousticanalysisdata['peak_t0'][evn, 0]
            p2_pk = acousticanalysisdata['peak_t0'][evn, 1]
            dpk.append(abs(p1_pk - p2_pk))
            rise_times_1.append(p1_pk-float(split_line[at00ind]))
            rise_times_2.append(p2_pk-float(split_line[at01ind]))
    gc.collect()
    i+=1
   
plt.figure()
plt.scatter(Vbaseline2, delta_at0)
plt.xlabel('piezo 2 baseline')
plt.ylabel(r'$|\Delta t_0|$ [s]')
plt.show()

plt.figure()
plt.scatter(Vbaseline1, delta_at0)
plt.xlabel('piezo 1 baseline')
plt.ylabel(r'$|\Delta t_0|$ [s]')
plt.show()
    
plt.figure()
plt.hist(rise_times_1, int(np.ceil(np.sqrt(len(rise_times_1)))))
plt.xlabel('piezo 1 rise time')
plt.show()

plt.figure()
plt.hist(rise_times_2, int(np.ceil(np.sqrt(len(rise_times_2)))))
plt.xlabel('piezo 2 rise time')
plt.show()

    
plt.figure()
plt.scatter(triggerwindowenergyp2, delta_at0)
plt.title('only bubble events')
plt.xlabel('piezo 2 energy in trigger window', fontsize=18)
plt.ylabel(r'$|\Delta t_0|$ [s]', fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-6, 1])

plt.show()

plt.figure()
plt.scatter(dpk, delta_at0)
plt.title('only bubble events')
plt.xlabel(r'$|\Delta t_{rise}|$ [s]', fontsize=18)
plt.ylabel(r'$|\Delta t_0|$ [s]', fontsize=18)

plt.show()