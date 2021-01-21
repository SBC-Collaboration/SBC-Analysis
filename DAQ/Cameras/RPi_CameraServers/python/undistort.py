import numpy as np
import cv2
import os

os.chdir("./Captures")
k = input("Enter the name of file to you want to undistort: ")
img = cv2.imread(k)
    

# Define Camera Matrix
mtx = np.loadtxt('cameraMatrix.txt', delimiter=',')


# Define distortion coefficients
dist = np.loadtxt('cameraDistortion.txt', delimiter=',')

# Undistort
dst = cv2.undistort(img, mtx, dist, None, mtx)
cv2.imwrite(k + 'Undistorted.png', dst)
print(f"Image {k} successfully undistorted")