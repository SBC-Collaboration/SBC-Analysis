#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:11:33 2020

@author: pi
"""

import numpy as np

results = np.subtract(background,current)
counter1 = np.count_nonzero(results>adc_threshold)
counter2 = np.count_nonzero(results<-1*adc_threshold)
counter = counter1 + counter2
if(counter>pix_threshold):
    feature_detect = True
if(feature_detect):
    for i in range(max_frames):
        im = Image.fromarray(ls[i])
        im = im.convert("L")
        im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
        print("images saved")
  