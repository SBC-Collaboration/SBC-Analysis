from capture import capture
import arducam_mipicamera as arducam
import v4l2
from datetime import datetime
import numpy as np
from PIL import Image

def set_controls(camera):
    try:
        camera.software_auto_exposure(enable=True)
    except Exception as e:
        print(e)
        
if __name__ == "__main__":
    try:
        camera = arducam.mipi_camera()
        camera.init_camera()
        print("camera open")
        camera.set_resolution(1920,1080)
        print("res set")
        camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
        camera.set_control(v4l2.V4L2_CID_HFLIP,1)
        set_controls(camera)
        ls = [None]*100
        for i in range(100):
            ls[i] = np.zeros((800,1280))
        ringBuf = tuple(ls)
    
        i = 0
        cont= True
        while(cont):
            
            try:
                if(i==100):
                    break
                else:
                    frame = camera.capture(encoding="raw")
                    buff1 = ringBuf[i]
                    buff2= frame.as_array.reshape(800,1280)
                    print("capture")
                    np.copyto(buff1,buff2)
                i = i+1
            except KeyboardInterrupt:
                break
        camera.close_camera()
        print("camera close")
        im = Image.fromarray(ringBuf[1])
        im = im.convert("L")
        im.save("test.png")
        for i in range(100):
            im = Image.fromarray(ringBuf[i])
            im = im.convert("L")
            im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
        print("images saved")
    except KeyboardInterrupt:
        print("ending image capture")