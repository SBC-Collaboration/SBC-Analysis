#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 02:49:57 2020

@author: David Velasco
"""
import numpy as np
from PIL import Image
#import RPi.GPIO as GPIO



ls = [None]*6
for i in range(6):
    ls[i] = np.zeros((1280,800))


for i in range(6):
    ls[i] = np.asarray(Image.open("/Users/dvelasco/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")) 
    

motion_trigger = 18
#GPIO.setmode(GPIO.BOARD)
#Need to confirm this number.....
#GPIO.setup(motion_trigger,GPIO.OUT,initial=GPIO.LOW)      
adc_threshold = 200
pix_threshold = 1000
background = ls[1]
current = ls[1]
for i in range(6):
    print(i)
    if(i==0):
        background = ls[i]
    else:
        current = ls[i]
        abs_diff = np.abs(np.subtract(current,background)) 
        counter = 0 
        for row in abs_diff:
            for pixel in row:
                if(pixel>adc_threshold):
                    counter +=1
        if(counter>pix_threshold):
              im = Image.fromarray(current) 
              im.show()
              #GPIO.output(motion_trigger, GPIO.HIGH)              
    i = i+1