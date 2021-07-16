#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 15:16:49 2021

@author: bressler
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 07:47:36 2019

@author: bressler
"""


import SBCcode as sbc
import numpy as np
import scipy
import matplotlib.pyplot as plt
from collections import Counter
import runlistscatalogue as rlc
from PICOcode.REFPROP.SeitzModel import SeitzModel
from os import listdir
from os.path import isfile,join


def xeSoundSpeed(runs,mergedfilename,Qbins,FeedbackTransducer):
    if FeedbackTransducer == "PT4":
        pt = 3
    elif FeedbackTransducer == "PT6":
        pt = 5

    #with open("/nashome/b/bressler/sbcoutput/%s_merged_oldPMT.txt"%mergedfilename,"r") as fin:
    with open("/nashome/b/bressler/sbcoutput/%s_merged.txt"%mergedfilename,"r") as fin:
        data = fin.readlines()
        
    headers = data[0].split()
    runind = headers.index("run")
    eventind = headers.index("event")
    xind = headers.index("x")
    yind = headers.index("y")
    zind = headers.index("z")
    lagind = headers.index("lag")
    spectind = headers.index("PMTphe")
    blockedind = headers.index('isBlocked')
    nbubind = headers.index('nbub')
    x=[]
    y=[]
    z=[]
    evid = []
    lag = []
    spect = []
    isblocked = []
    nbub = []
    i=0
    for line in data:
        if i>0:
            split_line = line.split()
            evid.append(split_line[runind]+"-"+split_line[eventind])
            isblocked.append(float(split_line[blockedind]))
            x.append(float(split_line[xind]))
            y.append(float(split_line[yind]))
            z.append(float(split_line[zind]))
            lag.append(float(split_line[lagind]))
            spect.append(float(split_line[spectind]))
            nbub.append(float(split_line[nbubind]))
    
        i+=1

    counts = []
    nbubs = []
    didntpass = 0
    bubInfo = []
    eventcount = np.zeros(87)
    spectra = [[] for i in range(87)]

    LT = np.zeros(87)
    LT_by_Q = np.zeros(len(Qbins)-1)
    expand_times = [[] for i in range(87)]
    totbub = 0
    onebubcount = np.zeros(87)
    twobubcount = np.zeros(87)
    threebubcount = np.zeros(87)
    setpoints = []
    elt = 0
    temps = []
    
    passing_lag = []
    passing_z = []
    
    for run in runs:
        print(run)
        tcut = 0 # was 20
        z_low = -4
        z_high = 0
        preq = 1
        runrawpath = '/bluearc/storage/SBC-17-data/'+run
        events = [evnt for evnt in listdir(runrawpath) if not isfile(join(runrawpath,evnt))]
        Nevents = len(events)
        runreconpath = "/pnfs/coupp/persistent/grid_output/SBC-17/output/%s/"%run
        historyfilename = runreconpath+"HistoryAnalysis_%s.bin"%run
        history = sbc.DataHandling.ReadBinary.ReadBlock(historyfilename)
        getbubfile = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/HumanGetBub_%s.bin"%run
        c = sbc.DataHandling.ReadBinary.ReadBlock(getbubfile)

        eventn = c["ev"]
        count = Counter(eventn)
        edges = history["PressureEdge"][pt]
        centersp = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        #centers = Qvals
        e0 = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,0,"slowDAQ")
        T=np.mean(e0["slowDAQ"]["T1"])
        temps.append(T)
        sm = SeitzModel(list(edges),T,'xenon')
        edges_by_Q = sm.Q
        centers = [(edges_by_Q[i] + edges_by_Q[i+1])/2 for i in range(len(edges_by_Q)-1)]

        print("T=%f C"%T)
        
        for c in count:
            N = count[c]
            counts.append(N)
            """
            if N%2 != 0 and N != 1:
                print("lines isn't even")
                print(c)
                print(N)
            elif N%2 == 0:
                n = N/2
                nbubs.append(n)
                """
        for eventn in range(Nevents):
            try:
                indices_back = 30
                event_ID = run+"-"+str(eventn)
                mergedind = evid.index(event_ID)
                bl = isblocked[mergedind]
                light = spect[mergedind]
                ev_z = z[mergedind]
                if np.isnan(light):
                    light = 0
                ev_z = z[mergedind]
                ev_lag = lag[mergedind]
                #onebub = not np.isnan(runposreco["z"][0][int(eventn)])
                n_from_merged = nbub[mergedind]
                n = 0
                N = count[eventn]
                if N%2 == 0:
                    n = N/2
                    nbubs.append(n)
                
                if n_from_merged != n:
                    print('n and merged n dont match')
                    print("n here: %d"%n)
                    print("n from merged: %d"%n_from_merged)
                e = sbc.DataHandling.GetSBCEvent.GetEvent(runrawpath,eventn,"slowDAQ","event")
                t = e["slowDAQ"]["elapsed_time"]
                tcenters = [(t[i+1] + t[i])/2 for i in range(len(t)-1)]
                trig_time = t[list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
                pslope = np.diff(e["slowDAQ"][FeedbackTransducer])
                expstartind = list(pslope).index(min(list(pslope)))
                exp_start = tcenters[expstartind]
                trig_pressure = e["slowDAQ"][FeedbackTransducer][list(e["slowDAQ"]["TriggerOut"]).index(1.0)-indices_back]
                T = np.mean(e["slowDAQ"]["T1"])
                
                #print(Q[0])
                pset = e["event"]["Pset"]
                ev_lt = e["event"]['livetime']
                elt += ev_lt
                if pset not in setpoints:
                    setpoints.append(pset)
                    #print(pset)
                #print(len(history["PressureBins"][eventn]))
                """
                plt.figure()
                plt.plot([trig_time,trig_time],[0,200],label="Bubble time")
                plt.plot([exp_start,exp_start],[0,200],label="Expansion start")
                plt.plot([0,max(t)],[trig_pressure,trig_pressure],label="Bubble pressure")
                plt.plot(t,e["slowDAQ"]["PT6"],label="PT6")
                plt.plot(tcenters,pslope,label="Pressure slope")
                plt.xlabel("time",fontsize=18)
                plt.ylabel("signal",fontsize=18)
                plt.legend(fontsize=18)
                plt.show
                """
                if n == 1 and ev_lt > 20 and pset == 25:
                    totbub += 1
                    for i in range(len(centers)):
                        if trig_time > exp_start + tcut:
                            #print(trig_pressure)
                            #print(pset)
                            if (trig_pressure >= edges[i] 
                                and trig_pressure <= edges[i+1]) and abs(trig_pressure-pset)<preq/2:
                                onebubcount[i] += 1
                                eventcount[i] += 1
                                if not bl and not np.isnan(light):
                                    if light>0:
                                        spectra[i].append(light)
                                        passing_lag.append(ev_lag)
                                        passing_z.append(ev_z)
            except Exception as x:
                print(x)
                                                        
    plt.figure()
    plt.scatter(passing_z, passing_lag, 20)
    plt.xlabel('z [cm]', fontsize=19)
    plt.ylabel('lag [microseconds]', fontsize=19)
    plt.show()
    
    def gaussian(x,mu,sigma,amplitude):
        return amplitude * np.exp(-((x - mu) /(np.sqrt(2)* sigma))**2 )

    def lin(x,m,b):
        return m*x + b
    
    zbinedges = np.arange(-2.75,0.25,0.25)
    zpts = [(zbinedges[i]+zbinedges[i+1])/2 for i in range(len(zbinedges)-1)]
    zbinlags = [[] for i in range(len(zbinedges)-1)]
    print(zbinlags)
    for j in range(len(passing_lag)):
        for i in range(len(zbinedges)-1):
            if passing_z[j] >= zbinedges[i] and passing_z[j] <zbinedges[i+1] and passing_lag[j] > -275:
                zbinlags[i].append(passing_lag[j])
    """
    params, params_cov = scipy.optimize.curve_fit(fit_fun, ch1zs, ch1lags)
    
    fitx = np.arange(-3,0,0.1)
    ft = [fit_fun(fitx[i], params[0], params[1], params[2],
                  params[3], params[4]) for i in range(len(fitx))]
    
    ss = (-1/params[3])*10000
    """
    
    avgs = []
    avgs_std = []
    
    gausfitx = np.arange(-400,0,0.1)
    
    for j in range(len(zpts)):
        Nbins=int(np.ceil(np.sqrt(len(zbinlags[j]))))
        print(Nbins)
        if Nbins == 0:
            Nbins = 1
        if len(zbinlags[j])>0:
            print(zbinlags[j])
            #avgs.append(np.mean(zbinlags[j]))
            gparams = [np.mean(zbinlags[j]), np.std(zbinlags[j]), 10]
            """
            dat,bins = np.histogram(zbinlags[j],Nbins)
            bincenters = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
            try:
                gparams,gparams_cov = scipy.optimize.curve_fit(gaussian,bincenters,dat,
                                                                   p0=(-250,50,10))
            except Exception as ex:
                print(ex)
                gparams = [np.mean(zbinlags[j]), np.std(zbinlags[j]), 10]
                print("bin centered on z=%.2f uses numpy mean and std"%zpts[j])
                print(gparams)
            """
            avgs.append(gparams[0])
            avgs_std.append(gparams[1])
            
            gausfit = [gaussian(gausfitx[i],gparams[0],gparams[1],gparams[2]) for i in range(len(gausfitx))]
            
        
    lfitx = np.arange(-3,0,0.1)
    print(zpts)
    print(avgs)
    print(avgs_std)
    lparams, lparams_cov = scipy.optimize.curve_fit(lin, zpts[:-1], avgs[:-1], sigma=avgs_std[:-1], p0 = [50, -50])
    ss = (1/abs(lparams[0]))*10000
    ss_err = 10000*max(abs(1/abs(lparams[0]) - 1/(abs(lparams[0])+np.sqrt(lparams_cov[0][0]))), abs(1/abs(lparams[0]) - 1/(abs(lparams[0])-np.sqrt(lparams_cov[0][0])))) 
    plt.figure()
    plt.grid()
    plt.scatter(passing_z,passing_lag)
    plt.plot(lfitx,lparams[0]*lfitx + lparams[1], color="g",linewidth=4,
             label = r"fit line: $\Delta t_0 = (%.2f \pm %.2f)z +(%.2f \pm %.2f)$"%(lparams[0], np.sqrt(lparams_cov[0][0]),
                                        lparams[1], np.sqrt(lparams_cov[1][1])))
    #plt.plot(fitx,params[3]*fitx + params[4],color = 'g',linewidth=4,
    #         label = "Linear part of fit, dt = %fz+%f"%(params[3],params[4]))
    #plt.plot(fitx,ft,color='r',linewidth=5,
    #         label = "Fit: dt = %fcos(%fz+%f) + %fz +%f"%(params[0],
    #                                 params[1],params[2],params[3],params[4]))
    plt.errorbar(zpts, avgs, avgs_std, color='c', marker='o',markersize=10, linestyle='none')
    plt.xlabel("z position [cm]",fontsize=20)
    plt.ylabel(r"$\Delta t_0 = t_{PMT} - t_{Piezo}$ [$\mu$s]",fontsize=20)
    plt.title(r"Channel 0, $v_s=%.2f \pm %.2f$ m/s"%(ss,ss_err),fontsize=20)
    plt.legend(fontsize=18)
    plt.show
    
    plt.figure()
    plt.hist(spect,int(np.floor(np.sqrt(len(spect)))))
    plt.show


def main():
    runs = rlc.cfJuly6to11
    #runs=['20170628_9']
    xeSoundSpeed(runs,'cfJuly6to11',[0,1,1.5,2,2.5,3],"PT6")
if __name__=="__main__":
    main()