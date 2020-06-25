import arducam_mipicamera as arducam
import v4l2
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

def set_controls(camera):
    try:
        camera.software_auto_exposure(enable=True)
    except Exception as e:
        print(e)        
if __name__ == "__main__":
    try:
#        GPIO.setmode(GPIO.BOARD)
#        GPIO.setup(36,GPIO.OUT,initial=GPIO.LOW)
        camera = arducam.mipi_camera()
        camera.init_camera()
        print("camera open")
        camera.set_resolution(1280,800)
        print("res set")
        camera.set_mode(5)
#        for i in range(7):
#            camera.write_sensor_reg(regs[i][0],regs[i][1])
        camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
        camera.set_control(v4l2.V4L2_CID_HFLIP,1)
        camera.set_control(v4l2.V4L2_CID_EXPOSURE,1)
#        set_controls(camera)
        ls = [None]*100
        for i in range(100):
            ls[i] = np.zeros((800,1280))
#        ringBuf = tuple(ls)
        i = 0
        t_end = time.time()+1
        while(True):
            try:
                if(i==100):
                    i = 0
                else:
                    frame = camera.capture(encoding="raw")
#                    buff1 = ringBuf[i]
                    ls[i]= frame.as_array.reshape(800,1280)
                    print(ls[i])
                    print("capture" +str(i))
#                    np.copyto(buff1,buff2)
                    i = i+1
            except KeyboardInterrupt:
                break
        camera.close_camera()
        print("camera close")
        for i in range(100):
            im = Image.fromarray(ls[i])
            im = im.convert("L")
            im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
        print("images saved")
    except KeyboardInterrupt:
        print("ending image capture")