import arducam_mipicamera as arducam
import v4l2 #sudo pip install v4l2
import time



if __name__ == "__main__":
    try:
        camera = arducam.mipi_camera()
        print("Open camera...")
        camera.init_camera()
        print("Setting the resolution...")
        fmt = camera.set_resolution(1920, 1080)
        camera.software_auto_exposure(enable = True)
        camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
        camera.set_control(v4l2.V4L2_CID_HFLIP,1)
        print("Current resolution is {}".format(fmt))
        print("Start preview...")
        camera.start_preview(fullscreen = False, window = (0, 0, 1280, 720))
        cont = input("Press e to stop")
        if (cont == "e"):
            print("Stop preview...")
            camera.stop_preview()
            print("Close camera...")
            camera.close_camera()
    except Exception as e:
        print(e)
