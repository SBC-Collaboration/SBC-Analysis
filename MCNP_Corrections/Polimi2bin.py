#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:20:38 2024

@author: hep
"""
import numpy as np
from ..DataHandling.WriteBinary import WriteBinaryNtupleFile as wb

def Polimi2bin(inputfilename, outputfilename):
    data = np.loadtxt(inputfilename)
    d = dict(id=np.int32(data[:,0]),
             pos=np.float64(data[:,8:11]),
             Er=np.float64(data[:,6]*1.0e3),
             species=np.int32(data[:,4]*0+132.0)
             )
    wb(outputfilename, d)