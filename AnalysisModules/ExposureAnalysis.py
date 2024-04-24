#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 14:59:54 2024

@author: cdahl
"""

import numpy as np


def ExposureAnalysis(ev,
                     PressureEdge=np.linspace(0,88,89, dtype=np.float64),
                     TimeEdge=np.float64([0, 5, 10, 20, 30]),
                     BubbleTime=np.float64([-0.205, -0.15])
                     ):

    EventPressure=-1+np.zeros((9), dtype=np.float64)
    Exposure=-1+np.zeros((9,len(PressureEdge)-1,len(TimeEdge)+1), dtype=np.float64)

    default_output = dict(PressureEdge=PressureEdge,
                          TimeEdge=TimeEdge,
                          BubbleTime=BubbleTime,
                          Exposure=Exposure,
                          EventPressure=EventPressure
                          )
    try:
        if not (ev['event']['loaded'] and ev['slowDAQ']['loaded']):
            return default_output

        dt = np.mean(np.diff(ev['slowDAQ']['elapsed_time']))
        trig_ix = np.nonzero(np.diff(ev['slowDAQ']['TriggerLatch'])==1)[0]+1
        exp_ix = np.nonzero(np.diff(ev['slowDAQ']['TriggerLatch'])==-1)[-1]
        bubix = trig_ix - np.round(BubbleTime/dt)
        lt_end = ev['slowDAQ']['elapsed_time'][trig_ix]-0.5*dt
        lt_start = lt_end - ev['event']['Event_livetime']
        ev_start = ev['slowDAQ']['elapsed_time'][0]
        ev_end = ev['slowDAQ']['elapsed_time'][-1]

        time_histedges = np.concatenate((ev_start[None]-999,
                                         TimeEdge+lt_start,
                                         ev_end[None]+999))

        for i_p in range(9):
            P = ev['slowDAQ']['PT'+str(i_p+1)]
            EventPressure[i_p] = P[bubix]
            P[:exp_ix] = -99
            P[trig_ix:] = -99
            Exposure[i_p] = dt * np.histogram2d(P, ev['slowDAQ']['elapsed_time'],
                                                [PressureEdge, time_histedges])
    except:
        return default_output
    return