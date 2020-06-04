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

class PiezoTab(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master.notebook)
        self.master = master