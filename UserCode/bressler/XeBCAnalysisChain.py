#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 12:26:58 2019

@author: bressler
"""

import SBCcode as sbc
import runlistscatalogue as rlc
import SBCcode.UserCode.bressler.SBC_AP as ap
import SBCcode.UserCode.bressler.find_bubble_PMT_coincidence as coinc
import SBCcode.UserCode.bressler.reconmerger as merger
import XeBCBubMult as mlt


def run_chain(runs, datasetname, gain, feedbackPT, LEDfraction, acousticPMTwindow, handscannedbool):
    """
    datadir = '/bluearc/storage/SBC-17-data/'
    for run in runs:
        ap.XeBCAP(datadir,run,'/nashome/b/bressler/sbcoutput/',20e3,120e3)
    """    
    coinc.find_coincidence(runs, gain, handscannedbool, datasetname, acousticPMTwindow)
    
    merger.merge(runs, handscannedbool, datasetname, gain, acousticPMTwindow)
    
    
    mlt.mult(runs, datasetname,[0,1,1.5,2,2.5,3],feedbackPT, [(2**i)+0.5 for i in range(16)], 'powerof2bins', LEDfraction)

def main():
    runs = rlc.bgJune27and28
    datasetname = 'bgJune27and28'
    feedbackPT = "PT6"
    run_chain(runs, datasetname, 1e7, feedbackPT, 0.1289, [-600e-6, -150e-6], False)
    
if __name__=="__main__":
    main()