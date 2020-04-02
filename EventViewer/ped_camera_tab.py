import getpass
import linecache
import logging
import os
import time
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk

import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from CustomToolbar import CustomToolbar
import numpy as np
from PIL import Image, ImageChops, ImageOps, ImageTk, PngImagePlugin, ImageDraw
import scipy.signal
from SBCcode.DataHandling.GetSBCEvent import GetEvent
from SBCcode.DataHandling import ReadBinary
import re

DEFAULT_FONT=("Arial", 20)

class CameraTab(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master.notebook)
        self.master = master
        
        self.init_image_width = 400
        self.init_image_height = 625
        self.native_image_width = None
        self.native_image_height = None
        self.max_zoom = 3
        
        self.image_naming_conventions = ['cam0_image0.png', 'cam0image  0.bmp']
        self.image_orientations = ['0', '90', '180', '270']
        self.plc_temp_var = 'T3'
        self.images_relative_path = 'Images'
        self.image_naming_convention = self.image_naming_conventions[0]
        self.num_cams = 4
        self.image_orientation = self.image_orientations[0]
        self.first_frame = '30'
        self.init_frame = '50'
        self.last_frame = '70'
        
        self.canvases = []
        
    def create_widgets(self):
        for cam in range(0, self.master.num_cams):
            canvas = tk.Canvas(self.master.camera_tab)
#            canvas.bind('<ButtonPress-1>', self.on_button_press)
#            canvas.zoom = 0
#    
#            canvas.image = canvas.create_image(0, 0, anchor=tk.NW, image=None)
#            canvas.bottom_text = canvas.create_text(10, self.init_image_height - 25, anchor=tk.NW, text='', fill='red', font=DEFAULT_FONT)
            canvas.grid(row=0, column=1 * cam, columnspan=1, sticky='NW')
            canvas.cam = cam
            self.canvases.append(canvas)
