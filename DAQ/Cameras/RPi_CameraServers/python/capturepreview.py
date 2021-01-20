import v4l2
import arducam_mipicamera as arducam
from datetime import datetime


def capture(camera):
    frame = camera.capture(encoding = 'jpeg')
    d1 = datetime.now().strftime("%d-%m-%Y_%H:%M:%S")
    path = "/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+d1+".jpg"
    open(path, "wb")
    frame.as_array.tofile(path) 
    del frame  # removes frame from memory


try:
    camera = arducam.mipi_camera()
    print("Open camera...")
    camera.init_camera()
    print("Setting the resolution...")
    fmt = camera.set_resolution(1920, 1080)
    camera.software_auto_exposure(enable=True)
    camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
    camera.set_control(v4l2.V4L2_CID_HFLIP, 1)  
        # Turning off HFLIP (mirror) may make capture easier but will need to flip again after capture if fed to cameracalib.py
    print("Current resolution is {}".format(fmt))
    print("Start preview...")
    camera.start_preview(fullscreen=False, window=(0, 0, 1280, 720))
    usrinput = input("Press ENTER to capture or input 'e' to stop")
    while True:
        if usrinput == '':
            capture(camera)
            print('Image captured!')
            usrinput = input("Press ENTER to capture or input 'e' to stop")
        elif usrinput == 'e':
            print("Stop preview...")
            camera.stop_preview()
            print("Close camera...")
            camera.close_camera()
        else:
            print('Error: invalid input')
            usrinput = input("Press ENTER to capture or input 'e' to stop")
except Exception as e:
    print(e)