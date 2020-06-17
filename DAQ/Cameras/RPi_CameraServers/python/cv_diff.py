#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 20:03:02 2020

@author: dvelasco
"""
import numpy as np
from numpy import asarray
from PIL import Image

for i in range(2-1):
    # load the image
    image1 = Image.open("/Users/dvelasco/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png") 
    image2 = Image.open("/Users/dvelasco/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i+1)+".png") 
    # convert image to numpy array
    data1 = asarray(image1)
    data2 = asarray(image2)
    # summarize shape
    print(data1.shape)
    print(data2.shape)
    im = Image.fromarray(np.subtract(data2,data2))
    im.show()
