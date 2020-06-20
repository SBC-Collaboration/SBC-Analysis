import arducam_mipicamera as arducam
import v4l2 #sudo pip3 install v4l2
import time
import ctypes
from datetime import datetime
from capture import capture

def set_controls(camera):
    try:
        print("Enable Auto Exposure...")
        camera.software_auto_exposure(enable = True)
        time.sleep(2)
        
    except Exception as e:
        print(e)

def callback(data):
    buff = arducam.buffer(data)
    file = buff.userdata
    buff.as_array.tofile(file)
    return 0
if __name__ == "__main__":
    try:
        camera = arducam.mipi_camera()
        print("Open camera...")
        camera.init_camera()
        print("Setting the resolution...")
        fmt = camera.set_resolution(1920, 1080)
        print("Current resolution is {}".format(fmt))
        camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
        camera.set_control(v4l2.V4L2_CID_HFLIP,1)
        print("Start preview...")
        camera.start_preview(fullscreen = False, window = (0, 0, 1280, 720))
        preview = True
        set_controls(camera)
        camera.software_auto_exposure(enable=False)
        camera.set_control(v4l2.V4L2_CID_EXPOSURE,4)
        cont = True
        d1 = datetime.now().strftime("%m-%d-%Y_%H:%M:%S")
        path = "/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Video/"+d1+".h264"
        file = open(path, "wb")
        # Need keep py_object reference
        file_obj = ctypes.py_object(file)
        camera.set_video_callback(callback, file_obj)
        while cont:
            cmd = input()
            cmd = cmd.lower()
            if (cmd == "s"):
                camera.stop_preview()
                preview = False
            elif(cmd=="p" and not preview):
                camera.start_preview(fullscreen=False, window = (0,0,1280,720))
                preview = True
            elif (cmd == "c"):
                capture(camera)
            elif(cmd == "stop"):
                cont = False
        
        camera.set_video_callback(None, None)
        file.close()
        print("Stop preview...")
        camera.stop_preview()
        print("Close camera...")
        camera.close_camera()
    except Exception as e:
        print(e)
