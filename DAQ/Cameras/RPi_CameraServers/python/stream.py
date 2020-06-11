from capture import capture
import arducam_mipicamera as arducam
import v4l2
from datetime import datetime
import numpy as np
#import keyboard

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
    ringBuf = (np.zeros((1280,800)),) *100
    i = 0
    cont= True
    while(cont):
#        if keyboard.is_pressed("q"):
#            cont=False
        if(i==100):
            # i = 0
            break
        else:
            buff = ringBuf[i]
            buff = camera.capture(encoding= "raw").as_array.reshape(1280,800)
            print("capture")
        i = i+1
    camera.close_camera()
    print("camera close")
    for i in range(100):
        open("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".raw","wb")
        ringBuf[i].tofile("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".raw")
        
    
        

