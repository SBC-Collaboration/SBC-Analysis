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

import ped_camera_tab  

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

PngImagePlugin.MAX_TEXT_CHUNK = 2000
DEFAULT_FONT=("Arial", 20)
# verbosity = logging.DEBUG
verbosity = logging.INFO

class PopUpHandler(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
        self.setFormatter(logging.Formatter(fmt='%(message)s'))
    def emit(self, message):
        messagebox.showerror('error', self.format(message))

logger = logging.getLogger('ped')
formatter = logging.Formatter(fmt='%(asctime)s [%(levelname)s] %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
console = logging.StreamHandler()
console.setFormatter(formatter)

logger.addHandler(console)
logger.setLevel(verbosity)

pop_up = PopUpHandler()
pop_up.setLevel(logging.ERROR)
logger.addHandler(pop_up)

class Application(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.grid()
        
        ###for running on local machine, change these based on local file location to set the correct data directories and initial dataset
        self.raw_directory = '/bluearc/storage/SBC-17-data/'
        self.base_directory, end = re.compile('\\w*-\\w*-data').split(self.raw_directory)
        self.scan_directory = '/coupp/data/home/coupp/scan_output_SBC-17/'
        self.reco_directory = '/pnfs/coupp/persistent/grid_output/SBC-17/output/'
        self.ped_directory = os.getcwd()
        self.npy_directory = os.getcwd()
        self.config_file_directory = os.path.join(self.ped_directory, 'configs')
        
        self.num_cams = 2
        
        self.notebook = ttk.Notebook(master)
        self.camera_tab = ped_camera_tab.CameraTab(self)
        
        self.notebook.add(self.camera_tab, text='camera')
        self.notebook.grid(row=0, column=0, columnspan=5)
        
        self.camera_tab.create_widgets()

ROOT = tk.Tk()
ROOT.lift()
ROOT.attributes('-topmost', True)
ROOT.after_idle(ROOT.attributes, '-topmost', False)
WIDTH = 1300  # width for the Tk root
HEIGHT = 900  # height for the Tk root
X = 0  # x coordinate to place root window at
Y = 0  # y coordinate to place root window at
# set the dimensions of the screen and where it is placed
ROOT.geometry('%dx%d+%d+%d' % (WIDTH, HEIGHT, X, Y))
ROOT.title('PED')
# handle quitting main loop on window close
def on_closing():
    ROOT.quit()
APP = Application(ROOT)
ROOT.protocol('WM_DELETE_WINDOW', on_closing)
APP.mainloop()