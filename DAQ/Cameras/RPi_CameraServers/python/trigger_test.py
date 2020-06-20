from capture import capture
import arducam_mipicamera as arducam
import v4l2
from datetime import datetime
import numpy as np
from PIL import Image
import time
import RPi.GPIO as GPIO

regs = [[0x4F00, 0x01],
        [0x3030, 0x04],
        [0x303F, 0x01],
        [0x302C, 0x00],
        [0x302F, 0x7F],
        [0x3823, 0x30],
        [0x0100, 0x00],]

camera = arducam.mipi_camera()
camera.init_camera()
print("camera open")
camera.set_resolution(1920,1080)
print("res set")
    
for i in range(7):
    camera.write_sensor_reg(regs[i][0],regs[i][1])
        
camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
camera.set_control(v4l2.V4L2_CID_HFLIP,1)
print("Start preview...")
camera.start_preview(fullscreen = False, window = (0, 0, 1280, 720))
cont = input("Press e to stop")
if (cont == "e"):
    print("Stop preview...")
    camera.stop_preview()
    print("Close camera...")
    camera.close_camera()
camera.stop_preview()
camera.close_camera()    

    
        
    
