#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 14:39:35 2020

@author: bressler
"""

import SBCcode as sbc
import runlistscatalogue as rlc
import matplotlib.pyplot as plt
import numpy as np
from PICOcode.REFPROP.SeitzModel import SeitzModel
from collections import Counter
from array import array
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile
   
def main():

    Et = [0.9, 1.19, 1.89, 2.06] # keV

    counts = [[57, 586, 39, 317], [1, 33, 1, 19]]
    lt = [6338.901, 63398.915, 4273.055, 50433.389]
    n_multiplicities = 2
    k = [n_multiplicities for i in range(len(Et))] # multiplicities; only doing single and multiple, not splitting up multiples
    b = 5 # nuisance parameters: cf, Y-be, Bi-Be, bg with veto, bg without veto
    # need index saying what thing the parameter effects; normalization, background, threshold
    #each row of file should have 2x5x3 matrix maybe 3x5x2 ***ndarray*** in numpy
    exp = [0.002376, 0.02383, 0.001602, 0.01896] # simulated live time: 2.668e6 sec
    bg_rate = np.array([[0.000756, 0.0003755, 0.0006757, 0.0005279], [0, 2.35e-5, 1.35e-4, 0]], dtype=float)
    nuisance = [np.zeros([3,b,n_multiplicities]) for i in range(len(Et))]
    #   0th index: which threshold
    #   1th index: which kind of parameter (normalization, background, threshold)
    #   2th index: which parameter (Cf, Y-Be, Bi-Be, bg with veto, bg without veto)
    #   3th index: which multiplicity it affects
    for i in range(len(Et)):
        for j in range(n_multiplicities):
            nuisance[i][0,0,j] = 0.012 # Cf-252 data fractional uncertainty on source activity
            nuisance[i][0,1,j] = 0 # this will be for Y-Be
            nuisance[i][0,2,j] = 0 # this will be for Bi-Be
    nuisance[0][1,3,0] = 0.169 # background rate fractional uncertainty, 0.90 keV
    nuisance[0][1,3,1] = 0.169 # background rate fractional uncertainty, 0.90 keV
    
    nuisance[1][1,3,0] = 0.25 # background rate fractional uncertainty, 1.14 keV
    nuisance[1][1,3,1] = 0.25 # background rate fractional uncertainty, 1.14 keV
    
    nuisance[2][1,3,0] = 0.316 # background rate fractional uncertainty, 1.89 keV
    nuisance[2][1,3,0] = 0.316 # background rate fractional uncertainty, 1.89 keV
    
    nuisance[3][1,3,0] = 0.354 # background rate fractional uncertainty, 2.06 keV
    nuisance[3][1,3,1] = 0.354 # background rate fractional uncertainty, 2.06 keV
    
    dat = {"E_T": np.array(Et,dtype=float),
           "max_mult": np.array(k,dtype=int),
           "exp": np.array(exp,dtype=float),
           "lt": np.array(lt,dtype=float),
           "bkg_rate": np.transpose(np.array(bg_rate,dtype=float)),
           "counts": np.transpose(np.array(counts,dtype=int)),
           "nuisance": np.array(nuisance,dtype=float)}
    print(dat)
    keylist = list(dat.keys())
    print(dat[keylist[0]])
    print([dat[k].shape[0] == dat[keylist[0]].shape[0]
                   for k in keylist])
    WriteBinaryNtupleFile("/nashome/b/bressler/sbcoutput/cfNRdata.bin",dat)
    #with open("/nashome/b/bressler/sbcoutput/NRCfAugust4and5.bin","wb") as outfile:
        #need to learn to write this stuff
        #dat=bytearray([Et,k,exp,lt,bg_rate,counts,nuisance])
        #outfile.write(dat)
    
    
if __name__=="__main__":
    main()