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

#LED_INDUCED_DEAD_TIME_FRACTION= 0.1289


def mult(runs,mergedfilename,Qbins,FeedbackTransducer, lightBins, binsname, LED_INDUCED_DEAD_TIME_FRACTION):
    if FeedbackTransducer == "PT4":
        pt = 3
    elif FeedbackTransducer == "PT6":
        pt = 5
    spectrumFile = open('/nashome/b/bressler/sbcoutput/'+mergedfilename+"%s_coincident_spectrum.txt"%binsname,'w')
    separator = '\t'
    spectrumFile.write("lightBin\t"+separator.join(["Q=%.1f-%.1fkeV"%(Qbins[i],
                                                                      Qbins[i+1]) for i in range(len(Qbins)-1)])+"\n")
    spectrumFile_multiples = open('/nashome/b/bressler/sbcoutput/'+mergedfilename+"%s_coincident_spectrum_multibub.txt"%binsname,'w')
    spectrumFile_multiples.write("lightBin\t"+separator.join(["Q=%.1f-%.1fkeV"%(Qbins[i],
                                                                      Qbins[i+1]) for i in range(len(Qbins)-1)])+"\n")
    fractionFile = open('/nashome/b/bressler/sbcoutput/'+mergedfilename+"%s_coincident_spectrum_fraction.txt"%binsname,'w')
    fractionFile.write("lightBin\t"+separator.join(["Q=%.1f-%.1fkeV"%(Qbins[i],
                                                                      Qbins[i+1]) for i in range(len(Qbins)-1)])+"\n")
    
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
    twobub_spectra = [[] for i in range(87)]
    threebub_spectra = [[] for i in range(87)]
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
        
    for run in runs:
        print(run)
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
                SM = SeitzModel(trig_pressure,T,'xenon')
                Q = (SM.Q)
                #print(Q[0])
                pset = e["event"]["Pset"]
                ev_lt = e["event"]['livetime']
                elt += ev_lt
                if pset not in setpoints:
                    setpoints.append(pset)
                    #print(pset)
                #print(len(history["PressureBins"][eventn]))
                times = history["PressureBins"][eventn][pt][:]
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
                if n == 1 and ev_lt > 20:
                    totbub += 1
                    for i in range(len(centers)):
                        if trig_time > exp_start:
                            #print(trig_pressure)
                            #print(pset)
                            if (trig_pressure >= edges[i] 
                                and trig_pressure <= edges[i+1]) and abs(trig_pressure-pset)<preq/2:
                                onebubcount[i] += 1
                                eventcount[i] += 1
                                if not bl:
                                    spectra[i].append(light)
                                    #if light < 1:
                                        #print("%s-%d is passing single not blocked and light < 1"%(run, eventn))
                                bubInfo.append([n,Q,bl,light,ev_z])
                                #print("added %s-%d to bubInfo"%(run,eventn))
                                expand_times[i].append(ev_lt)
                            if abs(centersp[i]-pset)<preq/2:
                                LT[i] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)
                                for j in range(len(Qbins)-1):
                                    if Q >= Qbins[j] and Q <= Qbins[j+1]:
                                        LT_by_Q[j] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)
                    if Q>1.5 and Q<2 and abs(trig_pressure-pset)<preq/2 and trig_time > exp_start:
                        print(event_ID)
                        print(Q)
                    
                        
                        
                elif n == 2 and ev_lt > 20:
                    totbub += 1
                    for i in range(len(centers)):
                        if trig_time > exp_start:
                            if (trig_pressure >= edges[i] 
                                and trig_pressure <= edges[i+1]) and abs(trig_pressure-pset)<preq/2:
                                eventcount[i] += 1
                                twobubcount[i] += 1
                                expand_times[i].append(ev_lt)
                                bubInfo.append([n,Q,bl,light,ev_z])
                                if not bl:
                                    twobub_spectra[i].append(light)
                                print("%s-%d multibubble event: %d bubbles, %.2f phe"%(run, eventn, n, light))
                            if abs(centersp[i]-pset)<preq/2:
                                LT[i] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)
                                for j in range(len(Qbins)-1):
                                    if Q >= Qbins[j] and Q <= Qbins[j+1]:
                                        LT_by_Q[j] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)

                elif n >= 3 and ev_lt > 20:
                    for i in range(len(centers)):
                        if trig_time > exp_start:
                            if (trig_pressure >= edges[i] 
                                and trig_pressure <= edges[i+1]) and abs(trig_pressure-pset)<preq/2:
                                eventcount[i] += 1
                                threebubcount[i] += 1
                                expand_times[i].append(ev_lt)
                                bubInfo.append([n,Q,bl,light,ev_z])
                                if not bl:
                                    threebub_spectra[i].append(light)
                                print("%s-%d multibubble event: %d bubbles, %.2f phe"%(run, eventn, n, light))
                            if abs(centersp[i]-pset)<preq/2:
                                LT[i] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)
                                for j in range(len(Qbins)-1):
                                    if Q >= Qbins[j] and Q <= Qbins[j+1]:
                                        LT_by_Q[j] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)

                else:
                    if n>0:
                        didntpass += 1
                    for i in range(len(centers)):
                        if trig_time > exp_start and ev_lt > 20:
                            if (trig_pressure >= edges[i] 
                                and trig_pressure <= edges[i+1]) and abs(trig_pressure-pset)<preq/2:
                                eventcount[i] += 1
                                expand_times[i].append(ev_lt)
                                bubInfo.append([n,Q,bl,light,ev_z])
                            if abs(centersp[i]-pset)<preq/2:
                                LT[i] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)
                                for j in range(len(Qbins)-1):
                                    if Q >= Qbins[j] and Q <= Qbins[j+1]:
                                        LT_by_Q[j] += times[i]*(1-LED_INDUCED_DEAD_TIME_FRACTION)

                #if (ev_z > z_low and ev_z < z_high):                    
                    
        
            except Exception as x:
                print(x)
                #break
        
    blockedEvents = [x for x in bubInfo if x[2] == 1]
    print("fraction of LED-blocked bubbles: "+str(len(blockedEvents)/len([x for x in bubInfo if x[0]>0])))
    print("Number of single-bubble zero phe non-blocked events:" + str(len([x for x in bubInfo if x[0]==1 and x[2] == 0 and (x[3]>0 and x[3] <1)])))
    print("Number of single-bubble one phe non-blocked events:" + str(len([x for x in bubInfo if x[0]==1 and x[2] == 0 and (x[3]>1 and x[3] <2)])))
    print("Number of single-bubble two phe non-blocked events:" + str(len([x for x in bubInfo if x[0]==1 and x[2] == 0 and (x[3]>2 and x[3] <3)])))

    rateList1 = [onebubcount[i]/LT[i] for i in range(len(onebubcount))]
    rateErrList1h = [(0.5+np.sqrt(onebubcount[i]+0.25))/LT[i] for i in range(len(onebubcount))]
    rateErrList1l = [(-0.5+np.sqrt(onebubcount[i]+0.25))/LT[i] for i in range(len(onebubcount))]
    standard_error1 = [np.sqrt(onebubcount[i])/LT[i] for i in range(len(onebubcount))]


    rateList2 = [twobubcount[i]/LT[i] for i in range(len(twobubcount))]
    rateErrList2h = [(0.5+np.sqrt(twobubcount[i]+0.25))/LT[i] for i in range(len(twobubcount))]
    rateErrList2l = [(-0.5+np.sqrt(twobubcount[i]+0.25))/LT[i] for i in range(len(twobubcount))]
    
    rateList3 = [threebubcount[i]/LT[i] for i in range(len(threebubcount))]
    rateErrList3h = [(0.5+np.sqrt(threebubcount[i]+0.25))/LT[i] for i in range(len(threebubcount))]
    rateErrList3l = [(-0.5+np.sqrt(threebubcount[i]+0.25))/LT[i] for i in range(len(threebubcount))]
    
    rateListE = [eventcount[i]/LT[i] for i in range(len(eventcount))]
    rateErrListEh = [(0.5+np.sqrt(eventcount[i]+0.25))/LT[i] for i in range(len(eventcount))]
    rateErrListEl = [(-0.5+np.sqrt(eventcount[i]+0.25))/LT[i] for i in range(len(eventcount))]
    """
    plt.figure()
    q = [bub[1] for bub in bubInfo]
    plt.hist(q,int(np.ceil(np.sqrt(len(q)))))
    plt.xlabel('Seitz Threshold')
    plt.ylabel('Count')
    plt.show()
    
    #Threshold:
    plt.figure()
    plt.errorbar(centers,rateList1,[rateErrList1l,rateErrList1h],fmt='ro',label='Single Bubbles')
    plt.errorbar(centers,rateList2,[rateErrList2l,rateErrList2h],fmt='go',label='Double Bubbles')
    plt.errorbar(centers,rateList3,[rateErrList3l,rateErrList3h],fmt='bo',label='Triple Bubbles')
    plt.errorbar(centers,rateListE,[rateErrListEl,rateErrListEh],fmt='ko',label='All Events')
    #plt.yscale('log')
    plt.ylim([1e-4,0.1])
    plt.xlabel('Seitz Threshold [keV]',fontsize=18)
    plt.ylabel('Rate [Hz]', fontsize=18)
    plt.legend(fontsize=18)
    plt.grid()
    plt.show
    
    #single bubbles
    def line(x,m,b):
        return m*x + b
    
    QForFit = []
    RForFit = []
    errForFit = []
    for i in range(len(rateList1)):
        if (not np.isnan(rateList1[i])) and (not np.isinf(rateList1[i])) and rateList1[i] != 0:
            QForFit.append(centers[i])
            RForFit.append(rateList1[i])
            errForFit.append(standard_error1[i])

    popt,pcov = scipy.optimize.curve_fit(line,QForFit,RForFit,p0=[0.005,0.01],sigma=errForFit)
    plt.figure()
    plt.errorbar(QForFit,RForFit,yerr=errForFit,fmt="ro")
    xfine = np.linspace(min(centers),3,50)
    y1 = line(xfine,popt[0]+pcov[0,0]**0.5, popt[1]-pcov[1,1]**0.5)
    y2 = line(xfine,popt[0]-pcov[0,0]**0.5, popt[1]+pcov[1,1]**0.5)
    plt.plot(xfine,line(xfine,popt[0],popt[1]),'k-')
    plt.plot(xfine,y1,'k--')
    plt.plot(xfine,y2,'k--')
    plt.plot([0,5],[0,0],'b-')
    plt.fill_between(xfine,y1,y2,facecolor='gray',alpha=0.1)
    plt.ylabel('Rate [Hz]',fontsize=18)
    plt.xlabel('Seitz Threshold [keV]',fontsize=18)
    plt.title('Crosses x-axis at %f'%(-popt[1]/popt[0]))
    plt.show
    """
    
    totbub = len(nbubs)
    print("total number of bubbles: %d"%totbub)
    print("total number of events: %d" %len(counts))
    print("didn't pass: %d"%didntpass)
    m = max(nbubs)
    for i in range(len(nbubs)):
        if nbubs[i]>3:
            nbubs[i] = 3
            
    binedges = [0.5,1.5,2.5,3.5] #bins for bubble multiplicity
    bincenters = [1,2,3]
    plt.figure()
    ns, _ = np.histogram(nbubs,binedges)
    print("bubble counts: ")
    print(ns)
    plt.errorbar(bincenters,ns,np.sqrt(ns),fmt='o')
    plt.yscale('log')
    plt.title("Max Nbub: %d"%int(m),fontsize=20)
    plt.xlabel("Nbub",fontsize=18)
    plt.ylabel("counts",fontsize=18)
    plt.xlim(0.5,3.5)
    plt.grid()
    plt.show
    
    f = [ns[i]/totbub for i in range(len(ns))]
    ferrh = [(0.5+np.sqrt(ns[i]+0.25))/totbub for i in range(len(ns))]
    ferrl = [(-0.5+np.sqrt(ns[i]+0.25))/totbub for i in range(len(ns))]

    plt.figure()
    plt.errorbar(bincenters,f,[ferrl,ferrh],fmt='o')
    plt.xlabel("Nbub",fontsize=18)
    plt.ylabel("fraction of bubbles",fontsize=18)
    plt.xlim([0.5,3.5])
    plt.grid()
    plt.show
    
    """
    #Pressure
    plt.figure()
    plt.errorbar(centersp,rateList1,[rateErrList1l,rateErrList1h],fmt='ro',label='Single Bubbles')
    plt.errorbar(centersp,rateList2,[rateErrList2l,rateErrList2h],fmt='go',label='Double Bubbles')
    plt.errorbar(centersp,rateList3,[rateErrList3l,rateErrList3h],fmt='bo',label='Triple Bubbles')
    plt.errorbar(centersp,rateListE,[rateErrListEl,rateErrListEh],fmt='ko',label='All Events')
    plt.yscale('log')
    plt.xlabel('Pressure [psia]',fontsize=18)
    plt.ylabel('Rate [Hz]', fontsize=18)
    plt.legend(fontsize=18)
    plt.grid()
    plt.show
    
    plt.figure()
    hy = []
    hx = []
    for i in range(len(expand_times)):
        if len(expand_times[i])>0:
            for j in range(len(expand_times[i])):
                hx.append(centersp[i])
                hy.append(expand_times[i][j])

    plt.hist2d(hx,hy,bins=(len(centers),30),cmap='magma')
    plt.colorbar()
    plt.xlabel('PT6 [psia]',fontsize=18)
    plt.ylabel('Total time since expansion [s]',fontsize=18)
    plt.show
    
        
    hy = []
    hx = []
    for i in range(len(expand_times)):
        if len(spectra[i])>0:
            for j in range(len(spectra[i])):
                hx.append(centersp[i])
                hy.append(spectra[i][j])
    plt.figure()
    plt.hist2d(hx,hy,bins=(40,int(max(hy))),cmap='Greys')
    plt.colorbar()
    plt.xlabel('PT6 [psia]',fontsize=18)
    plt.ylabel('Light collected [phe]',fontsize=18)
    plt.yscale('symlog',linthreshy = 0.9)
    plt.show
    
    
    
    xedges = np.arange(88)
    #yedges = np.arange(int(max(hy))+1)
    bins = [10**i+0.5 for i in range(5)] # light bins
    #bins = np.arange(0.5,1+np.ceil(max(spect)))
    bins = np.insert(bins,0,0.5)
    bins=np.insert(bins,0,-0.5)
    bins=np.insert(bins,0,-1.5)
    H, xedges, yedges = np.histogram2d(hx,hy,bins=(xedges,bins))
    for i in range(len(xedges)-1):
        for j in range(len(H[i,:])):
            #print(LT[i])
            H[i,j] = H[i,j]/LT[i]
    print("max rate: %f"%np.amax(H))
            
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    X,Y = np.meshgrid(xedges,yedges)

    im=ax.pcolormesh(X,Y,H.T,cmap="magma")
    plt.yscale('symlog',linthreshy = 0.9)
    fig.colorbar(im,ax=ax)
    plt.show()
    """
    hy = []
    hx = []
    hy_multi = []
    hx_multi = []
    for bubble in bubInfo:
        if bubble[0]==1 and bubble[2] == 0:
        #if bubble[0]==1:
            hx.append(bubble[1]) # Seitz threshold
            hy.append(bubble[3]) # phe
        elif (bubble[0] == 2 or bubble[0] == 3):
            hx_multi.append(bubble[1])
            hy_multi.append(bubble[3])
    plt.figure()
    plt.hist2d(hx,hy,bins=(40,int(max(hy))),cmap='Greys')
    plt.colorbar()
    plt.xlabel('Q [keV]',fontsize=18)
    plt.ylabel('Light collected [phe]',fontsize=18)
    plt.yscale('symlog',linthreshy = 0.9)
    plt.show
    
    #print(temps)
    
    xedges = edges
    Qedges = np.zeros(np.shape(xedges))
    for i in range(1,len(xedges)):
        SM = SeitzModel(float(xedges[i]),np.mean(temps),'xenon')
        if SM is not None:
            Qedges[i]=SM.Q
    
    OArateQ = [(Qbins[i]+Qbins[i+1])/2 for i in range(len(Qbins)-1)]
    OArate_multi = []
    OArateErr_multi = []
    OArate_singles = []
    OArateErr_singles = []
    Qedges = Qbins    
    #print(Qbins)
    #yedges = np.arange(int(max(hy))+1)
    
    bins = lightBins # light bins
    #bins = np.arange(1,int(1+np.ceil(max(spect))))
    bins = np.insert(bins,0,0.5)
    bins=np.insert(bins,0,-0.5)
    #bins=np.insert(bins,0,-1.5)
    
    #bins = [-0.5, 2**9+0.5, 2**13+0.5]
    binc=[(bins[i+1]+bins[i])/2 for i in range(len(bins)-1)]
    binwidths = [bins[i+1]-bins[i] for i in range(len(bins)-1)]
    print(binwidths)
    H, Qedges, yedges = np.histogram2d(hx,hy,bins=(Qedges,bins)) # using hx and hy defined a few lines above in the previous plot section
    H_multi, Qedges, yedges_multi = np.histogram2d(hx_multi, hy_multi, bins=(Qedges, bins))
    multi_ns = H_multi.copy()
    ns = H.copy()
    errorbararray_upper = np.zeros_like(H)
    errorbararray_lower = np.zeros_like(H)
    errorbararray_multi_upper = np.zeros_like(H_multi)
    errorbararray_multi_lower = np.zeros_like(H_multi)
    print("total number of events in H: %d"%np.sum(H))
    print("total number of multiples in H_multi: %d"%np.sum(H_multi))
    for i in range(len(LT_by_Q)):
        OArate_singles.append(sum(H[i,:])/LT_by_Q[i])
        OArateErr_singles.append(np.sqrt(sum(H[i,:]))/LT_by_Q[i])
        for j in range(len(H[i,:])):
            H[i,j] = H[i,j]/LT_by_Q[i]
            if ns[i,j] > 10:
                errorbararray_upper[i,j]= np.sqrt(ns[i,j])/LT_by_Q[i]
                errorbararray_lower[i,j]= np.sqrt(ns[i,j])/LT_by_Q[i]
            else:
                errorbararray_upper[i,j]= (0.5+np.sqrt(ns[i,j]+0.25))/LT_by_Q[i]
                errorbararray_lower[i,j]= (-0.5+np.sqrt(ns[i,j]+0.25))/LT_by_Q[i]
            #print(H[i,j])
            #print(errorbararray[i,j])
            if np.isnan(H[i,j]):
                H[i,j]=0
                errorbararray_upper[i,j] = 0
                errorbararray_lower[i,j] = 0
        OArate_multi.append(sum(H_multi[i,:])/LT_by_Q[i])
        OArateErr_multi.append(np.sqrt(sum(H_multi[i,:]))/LT_by_Q[i])
        for j in range(len(H_multi[i,:])):
            H_multi[i,j] = H_multi[i,j]/LT_by_Q[i]
            if multi_ns[i,j] > 10:
               errorbararray_multi_upper[i,j]= np.sqrt(multi_ns[i,j])/LT_by_Q[i]
               errorbararray_multi_lower[i,j]= np.sqrt(multi_ns[i,j])/LT_by_Q[i]
            else:
               errorbararray_multi_upper[i,j]= (0.5+np.sqrt(multi_ns[i,j] + 0.25))/LT_by_Q[i]
               errorbararray_multi_lower[i,j]= (-0.5+np.sqrt(multi_ns[i,j] + 0.25))/LT_by_Q[i]
            #print(H[i,j])
            #print(errorbararray[i,j])
            if np.isnan(H_multi[i,j]):
                H_multi[i,j]=0
                errorbararray_multi_upper[i,j] = 0
                errorbararray_multi_lower[i,j] = 0
    #print(sum(H))
    #print(OArate)
    print("singles:")
    print(ns)
    print("multiples:")
    print(multi_ns)
    
    plt.figure()
    plt.errorbar(OArateQ,OArate_singles,OArateErr_singles,marker='v',linestyle='')
    plt.grid()
    plt.xlabel('Seitz Threshold [keV]')
    plt.ylabel('Rate [Hz]')
    plt.show()
    
    print("max rate: %f"%np.amax(H))
    print("number in onebubcount: %f"%np.sum(onebubcount))
    print(np.sum([rateList1[i] * LT[i] for i in range(len(rateList1)) if not np.isnan(rateList1[i])]))
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    ax.title.set_text('Rates [Hz/bin]')
    X,Y = np.meshgrid(Qedges,yedges)

    im=ax.pcolormesh(X,Y,H.T,cmap="Greys")
    plt.yscale('symlog',linthreshy = 0.9)
    plt.ylabel("Collected light [phe]")
    plt.xlabel("Seitz Threshold [keV]")
    fig.colorbar(im,ax=ax)
    
    ax2 = fig.add_subplot(1,2,2)
    ax2.title.set_text('Number of bubbles in each bin')
    im2=ax2.pcolormesh(X,Y,ns.T,cmap="Blues")
    plt.yscale('symlog',linthreshy = 0.9)
    plt.ylabel("Collected light [phe]")
    plt.xlabel("Seitz Threshold [keV]")
    fig.colorbar(im2,ax=ax2)
    #plt.show()
    
    fig = plt.figure()
    ax3 = fig.add_subplot(1,1,1)
    ax3.title.set_text("Rates per Scintillation: fraction")
    for k in range(len(Qedges)-1):
        print("%f-%f keV: %f bubbles in %f seconds"%(Qedges[k],Qedges[k+1],sum(ns[k,:]),LT_by_Q[k]))
        if LT_by_Q[k]>100:
            ax3.errorbar(binc,H[k,:]/OArate_singles[k],[(errorbararray_upper[k,:]/OArate_singles[k]),(errorbararray_lower[k,:]/OArate_singles[k])],marker='.',markersize=7,
                         linestyle='none',label="%f-%f keV"%(Qedges[k],Qedges[k+1]))
    for l in range(len(bins)-1):
        fractionFile.write("%.1f-%.1fphe\t"%(bins[l],bins[l+1]))
        for k in range(len(Qedges)-1):
            fractionFile.write("%f,%f,%f\t"%(H[k,l]/OArate_singles[k], errorbararray_upper[k,l]/OArate_singles[k], errorbararray_lower[k,l]/OArate_singles[k]))
        fractionFile.write('\n')
    plt.xscale('symlog',linthreshx=0.5)
    #plt.yscale('symlog',linthreshy=0.0001)
    plt.legend(fontsize=15)
    plt.xlabel("collected light [phe]")
    plt.ylabel("Fraction")
    #plt.ylim([0,0.1])
    plt.grid(which='both',axis='both')
    plt.show()
    
    fig = plt.figure()
    ax3 = fig.add_subplot(1,2,1)
    #ax3.title.set_text("Rates per Scintillation")
    ax4 = fig.add_subplot(1,2,2)
    for k in range(len(Qedges)-1):
        print("%f-%f keV: %f seconds"%(Qedges[k],Qedges[k+1],LT_by_Q[k]))
        if LT_by_Q[k]>100:
            ax3.errorbar(binc,H[k,:],[errorbararray_upper[k,:], errorbararray_lower[k,:]],marker='.',markersize=7,
                         linestyle='none',label="%f-%f keV"%(Qedges[k],Qedges[k+1]))
            Hnew = []
            errorbarnew_upper = []
            errorbarnew_lower = []
            for l in range(len(binwidths)):
                Hnew.append(np.divide(H[k,l],binwidths[l]))
                errorbarnew_upper.append(np.divide(errorbararray_upper[k,l],binwidths[l]))
                errorbarnew_lower.append(np.divide(errorbararray_lower[k,l],binwidths[l]))
            ax4.errorbar(binc, Hnew, [errorbarnew_upper, errorbarnew_lower],
                         marker='.', markersize=7, linestyle='none',
                         label='%f-%f keV'%(Qedges[k],Qedges[k+1]))
    for l in range(len(bins)-1):
        spectrumFile.write("%.1f-%.1fphe\t"%(bins[l],bins[l+1]))
        for k in range(len(Qedges)-1):
            spectrumFile.write("%f,%f,%f\t"%(H[k,l],errorbararray_upper[k,l], errorbararray_lower[k,l]))
        spectrumFile.write('\n')
    ax3.set_xscale('symlog',linthreshx=0.5)
    ax4.set_xscale('symlog',linthreshx=0.5)
    ax3.set_yscale('log')
    ax4.set_yscale('log')
    #plt.yscale('symlog',linthreshy=0.0001)
    ax3.legend(fontsize=15)
    ax4.legend(fontsize=15)
    plt.xlabel("collected light [phe]")
    ax3.set_ylabel("Rate per bin [Hz]")
    ax4.set_ylabel("rate per phe [Hz/phe]")
    #plt.ylim([0,0.1])
    ax3.grid(which='both',axis='both')
    ax4.grid(which='both',axis='both')
    plt.show()
    
    fig = plt.figure()
    ax3 = fig.add_subplot(1,2,1)
    ax3.title.set_text("Rates per Scintillation, multiples")
    ax4 = fig.add_subplot(1,2,2)
    for k in range(len(Qedges)-1):
        print("%f-%f keV: %f seconds"%(Qedges[k],Qedges[k+1],LT_by_Q[k]))
        if LT_by_Q[k]>100:
            ax3.errorbar(binc,H_multi[k,:],[errorbararray_multi_upper[k,:], errorbararray_multi_upper[k,:]],marker='.',markersize=7,
                         linestyle='none',label="%f-%f keV"%(Qedges[k],Qedges[k+1]))
            Hnew_multi = []
            errorbarnew_multi_upper = []
            errorbarnew_multi_lower = []
            for l in range(len(binwidths)):
                Hnew_multi.append(np.divide(H_multi[k,l],binwidths[l]))
                errorbarnew_multi_upper.append(np.divide(errorbararray_multi_upper[k,l],binwidths[l]))
                errorbarnew_multi_lower.append(np.divide(errorbararray_multi_lower[k,l],binwidths[l]))
            ax4.errorbar(binc, Hnew_multi, [errorbarnew_multi_upper, errorbarnew_multi_lower],
                         marker='.', markersize=7, linestyle='none',
                         label='%f-%f keV'%(Qedges[k],Qedges[k+1]))
    for l in range(len(bins)-1):
        spectrumFile_multiples.write("%.1f-%.1fphe\t"%(bins[l],bins[l+1]))
        for k in range(len(Qedges)-1):
            spectrumFile_multiples.write("%f,%f,%f\t"%(H_multi[k,l],errorbararray_multi_upper[k,l], errorbararray_multi_lower[k,l]))
        spectrumFile_multiples.write('\n')
    ax3.set_xscale('symlog',linthreshx=0.5)
    ax4.set_xscale('symlog',linthreshx=0.5)
    ax3.set_yscale('log')
    ax4.set_yscale('log')
    #plt.yscale('symlog',linthreshy=0.0001)
    ax3.legend(fontsize=15)
    ax4.legend(fontsize=15)
    plt.xlabel("collected light [phe]")
    ax3.set_ylabel("Rate per bin [Hz]")
    ax4.set_ylabel("rate per phe [Hz/phe]")
    #plt.ylim([0,0.1])
    ax3.grid(which='both',axis='both')
    ax4.grid(which='both',axis='both')
    plt.show()


    
    plt.figure()
    plt.scatter([x[1] for x in bubInfo],[x[3] for x in bubInfo])
    plt.yscale('symlog',linthreshy=0.5)
    plt.xlabel('Seitz Threshold [keV]')
    plt.ylabel('collected light [phe]')
    plt.show

    
    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel('Seitz Threshold (keV)',fontsize=18)
    ax1.set_ylabel('Live time (s)',fontsize=18)
    ax2.set_ylabel('Number of bubbles',fontsize=18)
    ax1.scatter(centers,LT,20,'r')
    ax2.scatter(centers,onebubcount,20,'b')
    plt.show
    
    

    
    spectrumFile.close() 
    spectrumFile_multiples.close()       
    fractionFile.close()

def main():
    runs = rlc.cfJuneminus50CCombined
    #runs=['20170628_9']
    mult(runs,'cfJuneminus50CCombined',[0,1,1.5,2,2.5,3],"PT6", [(2**i)+0.5 for i in range(14)], 'powerof2bins', 0.1289)
if __name__=="__main__":
    main()