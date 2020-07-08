#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:19:14 2020

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
import logging
from multiprocessing import Process,Queue

class Stage1(Process):
    """Acquire frames as fast as possible and send to next stage"""
    def __init__(self, oqueue,camera):
        super().__init__()
        # Pick up parameters and store in class variables
        self.oqueue = oqueue      # output queue
        self.camera = camera

    def run(self,):
        # Turn on logging
        logging.basicConfig(level=logging.DEBUG,
                        format='%(created).6f [%(levelname)s] Stage1 %(message)s',
                        filename='log-stage1.txt', filemode='w')
        logging.info('started')

        # Generate frames and send down pipeline
        for f in range(100):
            logging.debug('Generating frame %d',f)
            # Generate frame of random stuff
            frame = camera.capture(encoding="raw")
            ls[f] = np.ctypeslib.as_array(frame.buffer_ptr[0].data,shape=(800,1280))
            logging.debug('Forwarding frame %d',f)
            self.oqueue.put(frame)

class Stage2(Process):
    """Read frames from previous stage as fast as possible, process and send to next stage"""
    def __init__(self, iqueue):
        super().__init__()
        # Pick up parameters and store in class variables
        self.iqueue = iqueue      # input queue


    def run(self,):
        # Turn on logging
        logging.basicConfig(level=logging.DEBUG,
                        format='%(created).6f [%(levelname)s] Stage2 %(message)s',
                        filename='log-stage2.txt', filemode='w')
        logging.info('started')
        
        for f in range(100):
            # Wait for next frame
            frame = self.iqueue.get()
            logging.debug('Received frame %d', f)
            # Process frame ...
            np.subtract(ls[f-1],ls[f-2],out=results)
            counter = count_above(results)
            logging.debug('Forwarding frame %d', f)
            if(counter > pix_threshold):
                return True
                

def main():
    q1_2 = Queue(2)
    
    stages = []
    stages.append(Stage1(q1_2,camera))
    stages.append(Stage2(q1_2))
    
    for stage in stages:
        stage.start()
        
    for stage in stages:
        stage.join()
        
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
    ls = np.zeros((100,800,1280),dtype=np.uint8)
    results = np.zeros((800,1280),dtype=np.uint8)
    adc_threshold1 = np.uint8(3)
    pix_threshold = np.uint8(199)
    main()