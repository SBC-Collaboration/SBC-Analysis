from capture import capture
import arducam_mipicamera as arducam
import v4l2
from datetime import datetime
import numpy as np
from PIL import Image
import time
import RPi.GPIO as GPIO

GPIO.setmode(GPIO.BOARD)
GPIO.setup(36,GPIO.OUT,initial=GPIO.LOW)
while True:
    GPIO.output(36,1)
    time.sleep(0.00416)
    GPIO.output(36,0)
    time.sleep(0.00416)
     
    