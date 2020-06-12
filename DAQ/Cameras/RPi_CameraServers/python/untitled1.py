from capture import capture
import arducam_mipicamera as arducam
import v4l2
from datetime import datetime
import numpy as np
#import keyboard
from PIL import Image

def set_controls(camera):
    try:
        camera.software_auto_exposure(enable=True)
    except Exception as e:
        print(e)
        
if __name__ == "__main__":
    camera = arducam.mipi_camera()
    camera.init_camera()
    print("camera open")
    camera.set_resolution(1920,1080)
    print("res set")
    camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
    camera.set_control(v4l2.V4L2_CID_HFLIP,1)
    set_controls(camera)
    ring = (np.zeros(1024000,),)
    frame = camera.capture(encoding="raw")
    camera.close_camera()
    print(ring[0])
    buff = ring[0]
    buff1 = frame.as_array
    np.copyto(buff,buff1)
    print(buff)
    print(ring[0])
    
    open("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(1)+".raw","wb")
    ring[0].tofile("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(1)+".raw")
   
