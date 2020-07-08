#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:24:33 2020

@author: pi
"""

import arducam_mipicamera as arducam
import v4l2
import numpy as np
from PIL import Image
from time import sleep
#import RPi.GPIO as GPIO
from ctypes import *
import ctypes
from count import count_above
from multiprocessing import Process

i=0
ls = np.zeros((100,800,1280),dtype=np.uint8)
results = np.zeros((800,1280),dtype=np.uint8)
adc_threshold1 = np.uint8(3)
pix_threshold = np.uint8(199)
camera = arducam.mipi_camera()

counter =0

def capture():
    global ls
    global camera
    global i
    if(i==0):
        camera.init_camera()
        print("camera open")
        camera.set_resolution(1280,800)
        camera.set_mode(5)
        print("res set")
        #    for i in range(7):
        #        camera.write_sensor_reg(regs[i][0],regs[i][1])
        camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
        camera.set_control(v4l2.V4L2_CID_HFLIP,1)
        camera.set_control(v4l2.V4L2_CID_EXPOSURE,1)
        
    frame = camera.capture(encoding="raw")
    print(i)
    ls[i] = np.ctypeslib.as_array(frame.buffer_ptr[0].data,shape=(800,1280))
    
def proc():
    global ls
    global i
    global results
    global adc_threshold1
    global counter
    np.subtract(ls[i-1],ls[i-2],out=results)
    counter = count_above(results, adc_threshold1)
    
if __name__=="__main__":
    
    capture()
    i+=1
    while(i<4):
        if(i==100):
            i=0
        p1 = Process(target=capture)
        p1.start()
        p2 = Process(target=proc)
        p2.start()
        if counter>pix_threshold:
            break
        i+=1
    camera.close_camera()
