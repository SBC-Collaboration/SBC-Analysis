#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:07:49 2020

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
import mpipe

ls = np.zeros((100,800,1280),dtype=np.uint8)
results = np.zeros((800,1280),dtype=np.uint8)
adc_threshold1 = np.uint8(3)
pix_threshold = np.uint8(199)
camera = arducam.mipi_camera()
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

def cap(i,camera,ls):
    
    i,camera,ls=i,camera,ls
    frame = camera.capture(encoding="raw")
    ls[i] = np.ctypeslib.as_array(frame.buffer_ptr[0].data,shape=(800,1280))
    background = ls[i-2]
    current = ls[i-1]
    return background, current
    
def processing(input_data):
    background,current = input_data
    np.subtract(background, current,out=results)
    counter = count_above(results, adc_threshold1)
    return counter
    
    


if __name__ == "__main__":
    
    
    
    for f in range(2):
        frame = camera.capture(encoding="raw")
        ls[f] = np.ctypeslib.as_array(frame.buffer_ptr[0].data,shape=(800,1280))
        del frame
    i = 2
    cont = True
    while cont:
        stage1 = mpipe.OrderedStage(cap, 1)
        stage2 = mpipe.OrderedStage(processing, 1)
        pipe = mpipe.Pipeline(stage1.link(stage2))
        
        pipe.put(i,camera,ls)
        
        if(i==100):
            i=0
        else:
            i+=1
        for result in pipe.results():
            counter = result
            if counter >pix_threshold:
                cont = False
        
    for i in range(100):
        im = Image.fromarray(ls[i])
        im = im.convert("L")
        im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
    print("images saved")
    
        

        


   