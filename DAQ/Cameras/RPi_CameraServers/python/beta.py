#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 23:28:39 2020

@author: David Velasco and Shishir B.
"""
import arducam_mipicamera as arducam
import v4l2
import numpy as np
from PIL import Image
import RPi.GPIO as GPIO
from count import count_above
import time
#import ray
import threading
import json
#ray.init()

i=0
ls = np.zeros((100,800,1280),dtype=np.uint8)
results = np.zeros((800,1280),dtype=np.uint8)
adc_threshold1 = np.uint8(10)
pix_threshold = np.uint8(199)
camera = arducam.mipi_camera()
GPIO.setmode(GPIO.BOARD)
GPIO.setup(7,GPIO.IN)
GPIO.setup(13,GPIO.IN)
GPIO.setup(11,GPIO.OUT, initial=GPIO.LOW)
#GPIO.setup(10,GPIO.OUT,initial=GPIO.LOW)

start_proc = threading.Event()
next_frame = threading.Event()
start_cap = threading.Event()

counter =0
loop=0

def wait_sig():
    while(True):
        if(GPIO.input(7)):
            start_cap.set()
            break

#@ray.remote
def capture():
    global ls
    global camera
    global i
    global loop
    print("start1")
    if(i==0 and loop==0):
        camera.init_camera()
        print("camera open")
        camera.set_resolution(1280,800)
        #camera.set_mode(5)
        print("res set")
        #    for i in range(7):
        #        camera.write_sensor_reg(regs[i][0],regs[i][1])
        camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
        camera.set_control(v4l2.V4L2_CID_HFLIP,1)
        camera.set_control(v4l2.V4L2_CID_EXPOSURE,1)
    if(i<2 and loop==0):
        pass
    else:
        start_proc.set()
    frame = camera.capture(encoding="raw")
    ls[i] = np.ctypeslib.as_array(frame.buffer_ptr[0].data,shape=(800,1280))
    
    
#@ray.remote
def proc():
    global ls
    global i
    global results
    global adc_threshold1
    global counter
    while(True):
        start_proc.wait()
        np.subtract(ls[i-1],ls[i-2],out=results)
        counter = count_above(results, adc_threshold1)
        print("done",counter)
        start_proc.clear()
        next_frame.set()
        
if __name__=="__main__":        
    t_end = time.time()+1
    p2 = threading.Thread(target=proc).start()
    while(True):
        i=0
        print("Waiting to begin capture")
        threading.Thread(target=wait_sig).start()
        start_cap.wait()
        start_cap.clear()
        print("Loading settings")
        adc_threshold1 = 10
        pix_threshold = 199
        frames_after = 10
        print("Loaded settings")
 #       GPIO.output(10,1)
        capture()
        i+=1
        capture()
        i+=1
        while(True):
            try:
                #if(GPIO.input(13)):
                if(loop==1):
                    break
                if(i==100):
                    i=0
                    loop+=1
                
                capture()
                if counter>pix_threshold:
                    GPIO.output(11,1)
                    print("Detected motion on frame " + str(i))
                next_frame.wait()
                next_frame.clear()
                i+=1

            except KeyboardInterrupt:
                break
        
        for j in range(frames_after):
            i+=1
            if(i==100):
                i=0
            capture()
        camera.close_camera()      
        for i in range(100):
            im = Image.fromarray(ls[i])
            im = im.convert("L")
            im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
        print("images saved")
  #      GPIO.output(10,0)
