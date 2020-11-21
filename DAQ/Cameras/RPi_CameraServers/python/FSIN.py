#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 19:46:43 2020

@author: pi     
"""

#creates .jpg file in Captures folder in SBCcode/Cameras/RPi_CameraServers/python/Captures
#file is labelled with data and time when taken

import arducam_mipicamera as arducam
import v4l2
from datetime import datetime
import time

regs = [[0x4F00, 0x01],
        [0x3030, 0x04],
        [0x303F, 0x01],
        [0x302C, 0x00],
        [0x302F, 0x7F],
        [0x3823, 0x30],
        [0x0100, 0x00],]


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self, k):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        fps = k/elapsed_time
        self._start_time = None
        print(f"Elapsed Time: {elapsed_time:0.4f} seconds")
        print(f"Frames: {k:0.4f}")
        print(f"Frames per second: {fps:0.4f} fps")


def set_controls(camera):
    try:
        print("Enable Auto Exposure...")
        camera.software_auto_exposure(enable = True)
    except Exception as e:
        print(e)
def capture(camera):
    frame = camera.capture(encoding = 'jpeg')
    d1 = datetime.now().strftime("%d-%m-%Y_%H:%M:%S")
    path = "/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+d1+".jpg"
    open(path, "wb")
    frame.as_array.tofile(path)
            #Remove frame from memory
    del frame
    
            
if __name__ == "__main__":
        try:
            camera = arducam.mipi_camera()
            print("Open camera...")
            camera.init_camera()
            print("Setting the resolution...")
            fmt = camera.set_resolution(1920, 1080)
            camera.set_mode(5)
            camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
            camera.set_control(v4l2.V4L2_CID_HFLIP,1)
            camera.set_control(v4l2.V4L2_CID_EXPOSURE,1)
            print("Current resolution is {}".format(fmt))
            t = Timer()
            k = 0
            for i in range(7):
                camera.write_sensor_reg(regs[i][0],regs[i][1])
            try:
                t.start()
                while(True):
                    capture(camera)
                    k = k+1
            except KeyboardInterrupt:
                print('Interrupted')      
            t.stop(k)
            print("Close camera...")
            camera.close_camera()
            
        except Exception as e:
            print(e)
            