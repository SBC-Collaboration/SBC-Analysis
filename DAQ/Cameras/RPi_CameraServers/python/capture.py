import arducam_mipicamera as arducam
#import v4l2
import pypng 
import time
from datetime import datetime

def set_controls(camera):
    try:
        print("Enable Auto Exposure...")
        camera.software_auto_exposure(enable = True)
    except Exception as e:
        print(e)

if __name__ == "__main__":
    try:
        camera = arducam.mipi_camera()
        print("Open camera...")
        camera.init_camera()
        print("Setting the resolution...")
        fmt = camera.set_resolution(1920, 1080)
        print("Current resolution is {}".format(fmt))
        set_controls(camera)
        time.sleep(1)
        frame = camera.capture(encoding = 'jpeg')
        today = datetime.now()
        d1 = today.strftime("%d/%m/%Y_%H:%M:%S")
        image = png.from_array(frame.as_array())
        image.save(d1+".png")
       

        # Release memory
        del frame
        print("Close camera...")
        camera.close_camera()
    except Exception as e:
        print(e)
