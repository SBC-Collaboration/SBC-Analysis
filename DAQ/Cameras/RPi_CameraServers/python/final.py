#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 10:25:56 2020

@author: pi
"""
import arducam_mipicamera as arducam
import v4l2
import numpy as np
from PIL import Image
import time
#import RPi.GPIO as GPIO
from ctypes import *
import ctypes
from count import count_above

regs = [[0x4F00, 0x01],
        [0x3030, 0x04],
        [0x303F, 0x01],
        [0x302C, 0x00],
        [0x302F, 0x7F],
        [0x3823, 0x30],
        [0x0100, 0x00],]


if __name__ == "__main__":
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
    adc_threshold1 = np.uint8(3)
    pix_threshold = 199 #15
    max_frames = 100
#    ls = [None]*100
#    for i in range(100):
#        ls[i] = np.zeros((800,1280),dtype=np.uint8)
    ls = np.zeros((max_frames,800,1280),dtype=np.uint8)
    results = np.zeros((800,1280),dtype=np.uint8)
    background = np.zeros((800,1280),dtype=np.uint8)
    current = np.zeros((800,1280),dtype=np.uint8)# prefilled array
    i = 0
    feature_detect = False
    bu = (c_ubyte*1024000) * 100 
    bf = bu()
    loop=0
    t_end=time.time()+1
    frame= camera.capture(encoding="raw")
    frame.buffer_ptr[0].data=bf[0]
    while(time.time()<=t_end):
        try:
            if(i==100):
                i = 0
            else:
                t_start=time.time()
                frame = camera.capture(encoding="raw")
#                print(addressof(frame.buffer_ptr[0].data))
#                print(frame.buffer_ptr[0].data)
#                print(addressof(frame.buffer_ptr.contents))
                print(time.time()-t_start)
                t_start = time.time()
#                if(i==99):
#                    frame.buffer_ptr[0].data=bf[0]
#                else:
#                    frame.buffer_ptr[0].data=bf[i+1]
#                ls[i]=frame.as_array.reshape(800,1280)
                ls[i]=np.ctypeslib.as_array(frame.buffer_ptr[0].data,shape=(800,1280))
#                print(addressof(bf[i]),addressof(frame.buffer_ptr[0].data))
                
                del frame
                print(time.time()-t_start)
                print("capture" +str(i))
#                if(i==0):
#                    t_start=time.time()
#                    background = ls[0]
#                    print(time.time()-t_start)
#                else:
#                    t_start=time.time()
#                    current=ls[i]
#                    print(time.time()-t_start)
#                    t_start=time.time()
#                    np.subtract(background,current,out=results)
#                    print(time.time()-t_start)
#                    t_start=time.time()
#                    counter1 = count_above(results,adc_threshold1)
#                    print(time.time()-t_start)
##                    if(counter1>pix_threshold):
#                    t_start=time.time()
#                    background = current
#                    print(time.time()-t_start)
                i+=1
        except KeyboardInterrupt:
            break
    camera.close_camera()
    for i in range(100):
        im = Image.fromarray(ls[i])
        im = im.convert("L")
        im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
    print("images saved")




