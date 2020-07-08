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
import time
#import ray

#ray.init()

i=0
ls = np.zeros((100,800,1280),dtype=np.uint8)
results = np.zeros((800,1280),dtype=np.uint8)
adc_threshold1 = np.uint8(10)
pix_threshold = np.uint8(199)
camera = arducam.mipi_camera()
done1= False
done2=False

counter =0
loop=0

#@ray.remote
def capture():
    global ls
    global camera
    global i
    global done1
    global loop
    print("start1")
    if(i==0 and loop==0):
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
    done1= True    
    frame = camera.capture(encoding="raw")
    print(frame.as_array)
    print(i)
    ls[i] = np.ctypeslib.as_array(frame.buffer_ptr[0].data,shape=(800,1280))
    done1 = False
    
#@ray.remote
def proc():
    global ls
    global i
    global results
    global adc_threshold1
    global counter
    global done2
    global done1
    print("start2")
    while(True):
        if(done1):
            np.subtract(ls[i-1],ls[i-2],out=results)
            counter = count_above(results, adc_threshold1)
            done2 = True
            print("done",counter)
    
if __name__=="__main__":
    
    capture()
    i+=1
    capture()
    i+=1
    
    t_end = time.time()+1
    p2 = Process(target=proc)
    p2.start()
    while(time.time()<=t_end):
        try:
            if(i==100):
                i=0
                loop+=1
            capture()
#            t_start = time.time()
#            p2 = Process(target=proc)
#            p2.start()
#            capture()
#            while(not done1 and not done2):
#                1+2
#            done1 = False
#            done2=False
#    #        ray.get([capture.remote(),proc.remote()])
            if counter>pix_threshold:
                break
            i+=1
#            print(time.time()-t_start)
        except KeyboardInterrupt:
            break
    camera.close_camera()
    for i in range(100):
        im = Image.fromarray(ls[i])
        im = im.convert("L")
        im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
    print("images saved")
