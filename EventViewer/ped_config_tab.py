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

class ConfigTab(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master.notebook)
        self.master = master
        
        # initialize variables
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
        
        # setup frames within tab
        self.config_tab_main = tk.Frame(self, bd=5, relief=tk.SUNKEN)
        self.config_tab_main.grid(row=0, column=0, sticky='NW')

        self.config_tab_vars = tk.Frame(self, bd=5, relief=tk.SUNKEN)
        self.config_tab_vars.grid(row=4, column=0, sticky='NW')

        # frame 1
        self.dataset_label = tk.Label(self.config_tab_main, text='Dataset:')
        self.dataset_label.grid(row=0, column=0, sticky='WE')
        self.dataset_select = ttk.Combobox(self.config_tab_main, values=self.get_datasets())
        self.dataset_select.grid(row=0, column=1, sticky='WE')
        self.set_init_dataset()

        self.update_dataset_button = tk.Button(self.config_tab_main, text='Update Dataset', command=self.update_dataset)
        self.update_dataset_button.grid(row=0, column=2, sticky='NW')

        self.raw_directory_label = tk.Label(self.config_tab_main, text='Raw Directory:')
        self.raw_directory_label.grid(row=1, column=0, sticky='WE')
        self.raw_directory_entry = tk.Entry(self.config_tab_main, width=50)
        self.raw_directory_entry.insert(0, self.raw_directory)
        self.raw_directory_entry.grid(row=1, column=1, sticky='WE')

        self.scan_directory_label = tk.Label(self.config_tab_main, text='Scan Directory:')
        self.scan_directory_label.grid(row=2, column=0, sticky='WE')
        self.scan_directory_entry = tk.Entry(self.config_tab_main, width=12)
        self.scan_directory_entry.insert(0, self.scan_directory)
        self.scan_directory_entry.grid(row=2, column=1, sticky='WE')

        self.reco_directory_label = tk.Label(self.config_tab_main, text='Reco Directory:')
        self.reco_directory_label.grid(row=3, column=0, sticky='WE')
        self.reco_directory_entry = tk.Entry(self.config_tab_main, width=12)
        self.reco_directory_entry.insert(0, self.reco_directory)
        self.reco_directory_entry.grid(row=3, column=1, sticky='WE')

        self.reco_version_label = tk.Label(self.config_tab_main, text='Reco Version:')
        self.reco_version_label.grid(row=4, column=0, sticky='WE')
        self.reco_version_combobox = ttk.Combobox(self.config_tab_main, values=['current', 'devel'])
        self.reco_version_combobox.grid(row=4, column=1, sticky='WE')

        self.npy_directory_label = tk.Label(self.config_tab_main, text='NPY Directory:')
        self.npy_directory_label.grid(row=5, column=0, sticky='WE')
        self.npy_directory_entry = tk.Entry(self.config_tab_main, width=12)
        self.npy_directory_entry.insert(0, self.npy_directory)
        self.npy_directory_entry.grid(row=5, column=1, sticky='WE')

        self.update_directory_button = tk.Button(self.config_tab_main, text='Update Directories',
                                                 command=self.update_directories)
        self.update_directory_button.grid(row=1, column=2, sticky='WE')
        
        # add quit button
        self.quit_button = tk.Button(self.config_tab_main, text='Quit PED', command=ROOT.destroy)
        self.quit_button.grid(row=2, column=2, sticky='WE')
        
        # frame 2
        self.manual_ped_config_directory_label = tk.Label(self.config_tab_vars, text='Manual Path to ped_config.txt')
        self.manual_ped_config_directory_label.grid(row=0, column=0, sticky='WE')
        self.manual_ped_config_directory = tk.Entry(self.config_tab_vars, width=70)
        self.manual_ped_config_directory.grid(row=0, column=1, sticky='NW')

        self.update_vars_config_button = tk.Button(self.config_tab_vars, text='Update Config Directory',
                                                   command=self.new_config_update)
        self.update_vars_config_button.grid(row=0, column=2, sticky='NW')

        self.update_vars_config_button = tk.Button(self.config_tab_vars, text='Update Config Values',
                                                   command=self.update_vars_config)
        self.update_vars_config_button.grid(row=10, column=2, sticky='NW')

        self.ped_config_file_path = tk.Label(self.config_tab_vars, text='Path to ped_config.txt file:')
        self.ped_config_file_path.grid(row=1, column=0, sticky='WE')
        self.ped_config_file_path_combobox = ttk.Combobox(self.config_tab_vars, width=60)
        self.ped_config_file_path_combobox.grid(row=1, column=1, sticky='WE')

        self.update_config_combobox_button = tk.Button(self.config_tab_vars, text='Update Config File',
                                                       command=self.config_combobox_update)
        self.update_config_combobox_button.grid(row=1, column=2, sticky='NW')

        self.plc_temp_config_label = tk.Label(self.config_tab_vars, text='PLC temperature var:')
        self.plc_temp_config_label.grid(row=2, column=0, sticky='WE')
        self.plc_temp_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.plc_temp_config_entry.insert(0, self.plc_temp_var)
        self.plc_temp_config_entry.grid(row=2, column=1, sticky='WE')

        self.relative_path_to_images_config_label = tk.Label(self.config_tab_vars, text='Relative path to images:')
        self.relative_path_to_images_config_label.grid(row=3, column=0, sticky='WE')
        self.relative_path_to_images_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.relative_path_to_images_config_entry.insert(0, self.images_relative_path)
        self.relative_path_to_images_config_entry.grid(row=3, column=1, sticky='WE')

        self.image_naming_convention_label = tk.Label(self.config_tab_vars, text='Image naming convention:')
        self.image_naming_convention_label.grid(row=4, column=0, sticky='WE')
        self.image_naming_convention_select = ttk.Combobox(self.config_tab_vars, values=self.image_naming_conventions)
        self.image_naming_convention_select.insert(0, self.image_naming_convention)
        self.image_naming_convention_select.grid(row=4, column=1, sticky='WE')

        self.num_cams_config_label = tk.Label(self.config_tab_vars, text='Number of cameras:')
        self.num_cams_config_label.grid(row=5, column=0, sticky='WE')
        self.num_cams_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.num_cams_config_entry.insert(0, self.num_cams)
        self.num_cams_config_entry.grid(row=5, column=1, sticky='WE')

        self.image_orientation_config_label = tk.Label(self.config_tab_vars, text='Image orientation:')
        self.image_orientation_config_label.grid(row=6, column=0, sticky='WE')
        self.image_orientation_select = ttk.Combobox(self.config_tab_vars, values=self.image_orientations)
        self.image_orientation_select.insert(0, self.image_orientation)
        self.image_orientation_select.grid(row=6, column=1, sticky='WE')

        self.first_frame_config_label = tk.Label(self.config_tab_vars, text='first frame:')
        self.first_frame_config_label.grid(row=7, column=0, sticky='WE')
        self.first_frame_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.first_frame_config_entry.insert(0, self.first_frame)
        self.first_frame_config_entry.grid(row=7, column=1, sticky='WE')

        self.init_frame_config_label = tk.Label(self.config_tab_vars, text='trig frame:')
        self.init_frame_config_label.grid(row=8, column=0, sticky='WE')
        self.init_frame_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.init_frame_config_entry.insert(0, self.init_frame)
        self.init_frame_config_entry.grid(row=8, column=1, sticky='WE')

        self.last_frame_config_label = tk.Label(self.config_tab_vars, text='last frame:')
        self.last_frame_config_label.grid(row=9, column=0, sticky='WE')
        self.last_frame_config_entry = tk.Entry(self.config_tab_vars, width=30)
        self.last_frame_config_entry.insert(0, self.last_frame)
        self.last_frame_config_entry.grid(row=9, column=1, sticky='WE')