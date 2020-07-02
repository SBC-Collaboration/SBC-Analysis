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
import ctypes
from count import count_above

regs = [[0x4F00, 0x01],
        [0x3030, 0x04],
        [0x303F, 0x01],
        [0x302C, 0x00],
        [0x302F, 0x7F],
        [0x3823, 0x30],
        [0x0100, 0x00],]


if __name__=="__main__":
    camera = arducam.mipi_camera()
    camera.init_camera()
    
    print("camera open")
    camera.set_resolution(1280,800)
    camera.set_mode(5)
    print("res set")
    #for i in range(7):
        #   camera.write_sensor_reg(regs[i][0],regs[i][1])
    camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
    camera.set_control(v4l2.V4L2_CID_HFLIP,1)
    camera.set_control(v4l2.V4L2_CID_EXPOSURE,4)
    #adc_threshold = np.uint8(3)
    pix_threshold = 100000 #15
    max_frames = 100
    #ls = [None]*100
    #for i in range(100):
    #    ls[i] = np.zeros((800,1280),dtype=np.uint8)
    ls = np.zeros((max_frames,800,1280), dtype = np.uint8)
    results = np.zeros((800,1280),dtype=np.uint8)
    i = 0
    feature_detect = False
    loop=0
    while(True):
        try:
            if(i==100):
                i = 0
            else:
                frame = camera.capture(encoding="raw")
#                    buff1 = ringBuf[i]
                print("capture" +str(i))
#                    np.copyto(buff1,buff2)
                if(i==0):
                    ls[i] = frame.as_array.reshape(800,1280)
                    print(str(i))
                else:
                    ls[i] = frame.as_array.reshape(800,1280)
                    t_start = time.time()
                    print(str(i))
                    results = np.subtract(ls[i-1],ls[i])
                    print(time.time()-t_start)
                    t_start=time.time()
            #        counter = (results>adc_threshold1).sum()
                    t_start=time.time()    
                    adc_threshold = np.uint8(np.sqrt(np.mean(ls[i-1])))     
                    print(time.time()-t_start)
                    t_start=time.time()
                    counter1 = count_above(results,adc_threshold)
            #        counter1 = (results>adc_threshold1).sum()
            #        counter2 = (results<adc_threshold2).sum()
                    print(time.time()-t_start)
                    counter = counter1
                    
                    if(counter>pix_threshold):
                        break
                i+=1
        except KeyboardInterrupt:
            break
    camera.close_camera()
    for i in range(100):
        im = Image.fromarray(ls[i])
        im = im.convert("L")
        im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
        print("images saved")