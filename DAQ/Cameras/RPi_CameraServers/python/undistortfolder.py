import numpy as np
import cv2
import os

read_directory = '/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures'
os.chdir(read_directory)

# Define Camera Matrix
mtx = np.loadtxt('cameraMatrix.txt', delimiter=',')


# Define distortion coefficients
dist = np.loadtxt('cameraDistortion.txt', delimiter=',')


# make a new directory to save undistorted images to
try:
    os.mkdir('/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Undistorted_Images')
except FileExistsError:
    print('Undistorted_Images directory already exists')
else:
    print('Made new directory: Undistorted_Images')

# undistort each image    
for image in os.listdir(read_directory): 
    if image.endswith('.jpg'):
        img = cv2.imread(image)
        image_name = os.path.splitext(image)[0]
        dst = cv2.undistort(img, mtx, dist, None, mtx)
        cv2.imwrite('/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Undistorted_Images/'+ image_name + 'Undistorted.png', dst)
        print(f"{image_name} succesfully undistorted")


