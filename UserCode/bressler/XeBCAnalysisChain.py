#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 12:26:58 2019

@author: bressler
"""

import SBCcode as sbc
from runlistscatalogue import *
import SBCcode.UserCode.bressler.SBC_AP as ap
import SBCcode.UserCode.bressler.coincidentbubblescintillation as coinc
import SBCcode.UserCode.bressler.reconmerger as merger

runs = cfJune28to30

datadir = '/bluearc/storage/SBC-17-data/'

for run in runs:
    ap.XeBCAP(datadir,run,'/nashome/b/bressler/sbcoutput/',20e3,120e3)

coinc.zdependence(runs)

merger.merge(runs,'cfJune28to30')

