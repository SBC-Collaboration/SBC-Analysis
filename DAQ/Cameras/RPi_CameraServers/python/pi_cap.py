#Working code for capture.py to incorporate GPIO logic and work with Arduino echoing signals


#creates .jpg file in Captures folder in SBCcode/Cameras/RPi_CameraServers/python/Captures
#file is labelled with data and time when taken

import arducam_mipicamera as arducam
import v4l2
from datetime import datetime
import RPi.GPIO as GPIO
import threading

GPIO.setmode(GPIO.BOARD)
GPIO.setup(11,GPIO.IN) #Connected to A2 (A2 -> R11) for clock signal 



result = None


def wait_sig():
    if(GPIO.input(11)):
        global result 
        result = 42
        


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

def main():
    thread = threading.Thread(target=wait_sig)
    thread.start()
    thread.join()
    camera = arducam.mipi_camera()
    print("Open camera...")
    camera.init_camera()
    print("Setting the resolution...")
    fmt = camera.set_resolution(1920, 1080)
    camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
    camera.set_control(v4l2.V4L2_CID_HFLIP,1)
    print("Current resolution is {}".format(fmt))
    set_controls(camera)
#           sleep(1)
#            capture(camera)
    for i in range(100):
        capture(camera)
        print(i)
    print("Close camera...")
    camera.close_camera()
      
    
    

         
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
        
            