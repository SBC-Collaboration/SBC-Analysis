#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:46:43 2020

@author: bressler
"""

import SBCcode as sbc
import numpy as np
import runlistscatalogue as rlc
from dogammasmakebubbles import getRate

allxyzfname = "/pnfs/coupp/persistent/grid_output/SBC-17/output/SimpleXYZ_all.bin"
xyzf = sbc.DataHandling.ReadBinary.ReadBlock(allxyzfname)
print(xyzf.keys())
pressureRequirement = 2.5

print("----------------------------------------")
print("Test: Y-Be")
yberuns = rlc.YBeSept26to28
getRate(xyzf,yberuns,pressureRequirement)

print("---------------------------------")
print("background")
print("-45C, july")
bgm45runs = rlc.bgJuly11and12
getRate(xyzf,bgm45runs,pressureRequirement)
print("-45C, september")
bgm45septruns = rlc.bgSept27
getRate(xyzf,bgm45septruns,pressureRequirement)
print("-43C, october")
bgm43runs = rlc.bgOct2and3
getRate(xyzf,bgm43runs,pressureRequirement)


print("------------------------------------------")
print("Cobalt-57")
co57m45 = rlc.CoJuly5and6
getRate(xyzf,co57m45,pressureRequirement)

print("---------------------------------------------")
print("Y-Al")
print("Type 8850: -45")
yal8850m45runs = rlc.YAlSept26and27
getRate(xyzf,yal8850m45runs,pressureRequirement)
print("Type 8850: -43")
yal8850m43runs = rlc.YAlOct1
getRate(xyzf,yal8850m43runs,pressureRequirement)

print("-------------------------------------------")
print("BiAl")
bialruns = rlc.BiAlSept24and25
getRate(xyzf,bialruns,pressureRequirement)

print("-------------------------------------------")
print("Bi-207")
bim45runs  = rlc.BiJuly6
getRate(xyzf,bim45runs,pressureRequirement)


