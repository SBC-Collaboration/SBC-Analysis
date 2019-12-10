#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:08:39 2019

@author: bressler
"""

import SBCcode as sbc
from cobalt57spectrumattempt import spectrum
from coincidentbubblescintillation import zdependence

run = '20170803_2'
datapath = '/bluearc/storage/SBC-17-data/'
runpath = datapath+run

a,b,c = spectrum(runpath)

[goodz,diffs,coincspec,Ncoinc,ntotcoinc,totevents,totbub]=zdependence([run])