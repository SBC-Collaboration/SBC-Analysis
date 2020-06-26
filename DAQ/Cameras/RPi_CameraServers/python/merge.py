#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 

@author: pi
"""
import arducam_mipicamera as arducam
import v4l2
import numpy as np
from PIL import Image
import time
#import RPi.GPIO as GPIO
import ctypes


regs = [[0x4F00, 0x01],
        [0x3030, 0x04],
        [0x303F, 0x01],
        [0x302C, 0x00],
        [0x302F, 0x7F],
        [0x3823, 0x30],
        [0x0100, 0x00],]
 
#_____________________________________________________________________________  
class DynamicArray(object): 
    ''' 
    DYNAMIC ARRAY CLASS (Similar to Python List) 
    '''
      
    def __init__(self): 
        self.n = 0 # Count actual elements (Default is 0) 
        self.capacity = 1 # Default Capacity 
        self.A = self.make_array(self.capacity) 
          
    def __len__(self): 
        """ 
        Return number of elements sorted in array 
        """
        return self.n 
      
    def __getitem__(self, k): 
        """ 
        Return element at index k 
        """
        if not 0 <= k <self.n: 
            # Check it k index is in bounds of array 
            return IndexError('K is out of bounds !')  
          
        return self.A[k] # Retrieve from the array at index k 
          
    def append(self, ele): 
        """ 
        Add element to end of the array 
        """
        if self.n == self.capacity: 
            # Double capacity if not enough room 
            self._resize(2 * self.capacity)  
          
        self.A[self.n] = ele # Set self.n index to element 
        self.n += 1
    
    def insertAt(self,item,index): 
        """ 
         This function inserts the item at any specified index. 
        """
  
          
        if index<0 or index>self.n: 
            print("please enter appropriate index..") 
            return
          
        if self.n==self.capacity: 
            self._resize(2*self.capacity) 
              
          
        for i in range(self.n-1,index-1,-1): 
            self.A[i+1]=self.A[i] 
              
          
        self.A[index]=item 
        self.n+=1
  
  
          
    def delete(self): 
        """ 
        This function deletes item from the end of array 
        """
  
        if self.n==0: 
            print("Array is empty deletion not Possible") 
            return
          
        self.A[self.n-1]=0
        self.n-=1
          
          
          
      
    def removeAt(self,index): 
        """ 
        This function deletes item from a specified index.. 
        """        
  
        if self.n==0: 
            print("Array is empty deletion not Possible") 
            return
                  
        if index<0 or index>=self.n: 
            return IndexError("Index out of bound....deletion not possible")         
          
        if index==self.n-1: 
            self.A[index]=0
            self.n-=1
            return        
          
        for i in range(index,self.n-1): 
            self.A[i]=self.A[i+1]             
              
          
        self.A[self.n-1]=0
        self.n-=1
  
          
    def _resize(self, new_cap): 
        """ 
        Resize internal array to capacity new_cap 
        """
          
        B = self.make_array(new_cap) # New bigger array 
          
        for k in range(self.n): # Reference all existing values 
            B[k] = self.A[k] 
              
        self.A = B # Call A the new bigger array 
        self.capacity = new_cap # Reset the capacity 
          
    def make_array(self, new_cap): 
        """ 
        Returns a new array with new_cap capacity 
        """
        return (new_cap * ctypes.py_object)() 

#______________________________________________________________________________
def set_controls(camera):
    try:
        camera.software_auto_exposure(enable=True)
    except Exception as e:
        print(e)        
#______________________________________________________________________________
if __name__ == "__main__":
    try:
        #GPIO.setmode(GPIO.BOARD)
        #GPIO.setup(36,GPIO.OUT,initial=GPIO.LOW)
        camera = arducam.mipi_camera()
        camera.init_camera()
        camera.set_mode(5)
        print("camera open")
        camera.set_resolution(1280,800)
        print("res set")
        #for i in range(7):
         #   camera.write_sensor_reg(regs[i][0],regs[i][1])
        camera.set_control(v4l2.V4L2_CID_VFLIP, 1)
        camera.set_control(v4l2.V4L2_CID_HFLIP,1)
        camera.set_control(v4l2.V4L2_CID_EXPOSURE,4)
#        set_controls(camera)
        adc_threshold = 3
        pix_threshold = 1 #15
        max_frames = 120
        ls = np.zeros((max_frames,1280,800))
#      entries = range(1024000) # 1 million entries
        results = np.zeros((1280,800)) # prefilled array
        i = 0
        feature_detect = False 
        t_end = time.time()+1
        while(time.time()<t_end):
            try:
                if(i==max_frames):
                    for j in range(max_frames):
                        if(j==0):
                            background = ls[j]
                        else:
                            current = ls[j]
                            results = np.abs(np.subtract(background,current))
                            counter = np.count_nonzero(results>adc_threshold)
                            if(counter>pix_threshold):
                                feature_detect = True
                                break
                                #print("All Zeros") 
                            background = current
                        #print("capture" +str(j))
                    i= -1
                else:
                    frame = camera.capture(encoding="raw")
                    ls[i]=frame.as_array.reshape(1280,800)
                    print(i)
                    i +=1
            except KeyboardInterrupt:
                break
        camera.close_camera()
        if(feature_detect):
            for i in range(max_frames):
                im = Image.fromarray(ls[i])
                im = im.convert("L")
                im.save("/home/pi/SBCcode/DAQ/Cameras/RPi_CameraServers/python/Captures/"+str(i)+".png")
                print("images saved")
        print("camera close")
    except KeyboardInterrupt:
        print("ending image capture")
