#!/coupp/app/home/coupp/anaconda3/bin/python
'''pico-event-display

to run: python ped.py
may need to add to your paths:
 export PATH=/coupp/app/home/coupp/anaconda3/bin:$PATH
 export PYTHONPATH=/coupp/app/home/coupp/anaconda3/pkgs:$PYTHONPATH

v0.5.1 Karl August 2, 2017 -- Bug fixes and release on COUPP machine
v0.5.0 Karl July 7, 2017 -- Added dytran, config tab, bug fixes
v0.4.7 Orin Jan 10, 2017 -- Totally new way of handling merged_all using npy files from convert.py script, and some interface improvements
v0.4.6 Orin Dec 20, 2016 -- Now all piezo traces are available
v0.4.4 Orin Dec 20, 2016 -- Added fastDAQ_t0 to piezo trace and merged output vars, added buttons for quickly traversing the merged_all file by 1K lines at a time
v0.4.3 Orin Dec 19, 2016 -- Added reading in of fastDAQ data, drawing of piezo traces
v0.4.2 Orin Dec 18, 2016 -- fixed a bug when finding the crosshairs, optimized event traversal
v0.4.1 Orin Dec 17, 2016 -- changes to file handling to deal with changing merged_file, added crosshairs, added tab for piezo traces, put GUI elements in separate canvases
v0.3.8 Orin Nov 3, 2016 -- added cut2, and now stops searching for next event if fails to find event passing cut 1000 times in a row
v0.3.7 Orin Oct 27, 2016 -- add option to load merged_all.txt file, and to cut on a variable in the merged file
v0.3.5 Orin Oct 6, 2016 -- added 'first frame' and 'last frame' buttons, and made diff go to frame 30
v0.2.7 Orin Sep 4, 2016
v0.2.7-linted1: pylint compliance mods, PEP8 (w/120-char line limit) by Scott Sep 6, 2016
v1.0 Orin Sep 9, 2016 -- added handscanning and forward/back navigation through raw data folders

'''
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

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

PngImagePlugin.MAX_TEXT_CHUNK = 2000
# verbosity = logging.DEBUG
verbosity = logging.INFO

DEFAULT_FONT=("Arial", 20)

class PopUpHandler(logging.Handler):
    # x
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
        self.init_image_width = 400
        self.init_image_height = 625
        self.native_image_width = None
        self.native_image_height = None
        self.max_zoom = 3

        # raw_directory: where the data for an individual dataset is stored
        # base_directory: where the datasets are stored, created from raw_directory
        # scan_directory: where the handscans are stored
        # reco_directory: where the .npy files (raw_events.npy, reco_events.npy) and merged_all.txt
        # ped_directory: where the ped code is stored on the machine
        # config_file_directory: where the configuration files are placed, should be a folder in the same directory as ped_directory
        ### also change self.ped_config_file_path_var to append desired initial dataset rather than default

        # # for running on COUPP machine
        # self.raw_directory = '/bluearc/storage/30l-16-data/'
        # self.base_directory, end = re.compile('\\w*-\\w*-data').split(self.raw_directory)
        # self.scan_directory = '/coupp/data/home/coupp/scan_output_30l-16/'
        # self.reco_directory = '/coupp/data/home/coupp/recon/current/30l-16/output/'
        # self.ped_directory = '/coupp/data/home/coupp/PEDsvn/'
        # self.config_file_directory = os.path.join(self.ped_directory, 'configs')

        ###for running on local machine, change these based on local file location to set the correct data directories and initial dataset
        self.raw_directory = '/bluearc/storage/SBC-17-data/'
        self.base_directory, end = re.compile('\\w*-\\w*-data').split(self.raw_directory)
        self.scan_directory = '/coupp/data/home/coupp/scan_output_SBC-17/'
        self.reco_directory = '/pnfs/coupp/persistent/grid_output/SBC-17/output/'
        #self.reco_directory = '/coupp/data/home/coupp/recon/current/30l-16/output/'

        self.ped_directory = os.getcwd()
        self.npy_directory = os.getcwd()

        self.config_file_directory = os.path.join(self.ped_directory, 'configs')

        # #  for running on local machine, change these based on local file location to set the correct data directories and initial dataset
        # self.raw_directory = './30l-16-data'
        # self.base_directory, end = re.compile('\\w*-\\w*-data').split(self.raw_directory)
        # self.scan_directory = './scan_output_30l-16'
        # self.reco_directory = './30l-16-data/output'
        # self.ped_directory = "./"
        # self.config_file_directory = './configs'

        # Errors will be appended to this string and displayed for each event
        self.error = ''
        # Default config values
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
        self.piezo = 'Piezo7'
        self.ped_config_file_path_var = os.path.join(self.config_file_directory, 'SBC-17-ped_config.txt')
        self.load_config_values(self.ped_config_file_path_var)
        self.n_PMT_trig = tk.IntVar(value=0)
        self.current_PMT_event = None
        self.get_raw_events()
        self.source_button_var = tk.IntVar(value=-1)
        self.nbub_button_var = tk.IntVar(value=-1)
        self.do_handscan_checkbutton_var = tk.BooleanVar(value=False)
        self.draw_crosshairs_var = tk.BooleanVar(value=False)
        self.show_bubbles_var = tk.BooleanVar(value=True)
        self.invert_checkbutton_var = tk.BooleanVar(value=False)
        self.diff_checkbutton_var = tk.BooleanVar(value=False)
        self.antialias_checkbutton_var = tk.BooleanVar(value=True)
        self.load_dytran_checkbutton_var = tk.BooleanVar(value=False)
        self.piezo_plot_t0_checkbutton_var_top = tk.BooleanVar(value=False)
        self.piezo_plot_t0_checkbutton_var_bottom = tk.BooleanVar(value=False)
        self.piezo_plot_t0_checkbutton_vars = [self.piezo_plot_t0_checkbutton_var_top,
                                               self.piezo_plot_t0_checkbutton_var_bottom]
        self.dytran_plot_t0_checkbutton_var = tk.BooleanVar(value=False)
        self.load_fastDAQ_piezo_checkbutton_var_top = tk.BooleanVar(value=False)
        self.load_fastDAQ_piezo_checkbutton_var_bottom = tk.BooleanVar(value=False)
        self.load_fastDAQ_piezo_checkbutton_vars = [self.load_fastDAQ_piezo_checkbutton_var_top,
                                                    self.load_fastDAQ_piezo_checkbutton_var_bottom]
        self.isgoodtrigger_checkbutton_var = tk.BooleanVar(value=True)
        self.crosshairsgood_checkbutton_var = tk.BooleanVar(value=True)
        for i in range(9):
            self.grid_rowconfigure(i, weight=1)
        self.scanner_name = tk.StringVar()
        # PLC vars
        self.temp_label = tk.StringVar()
        # event.txt vars
        self.run_type = -1
        self.run_type_label = tk.StringVar()
        self.pset_label = tk.StringVar()
        self.te_label = tk.StringVar()
        # For the fastDAQ tab
        self.piezo_cutoff_low = 2000
        self.piezo_cutoff_high = 10000
        self.piezo_beginning_time = -.1
        self.piezo_ending_time = 0.0
        self.incremented_piezo_event = False
        self.piezo_timerange_checkbutton_var_top = tk.BooleanVar(value=False)
        self.piezo_timerange_checkbutton_var_bottom = tk.BooleanVar(value=False)
        self.piezo_timerange_checkbutton_vars = [self.piezo_timerange_checkbutton_var_top,
                                                 self.piezo_timerange_checkbutton_var_bottom]
        # For config window
        self.coupp_checkbutton_var = tk.BooleanVar(value=False)
        self.run = None
        self.event = None
        self.row_index = -1
        self.cuts = []
        self.selected_events = None
        self.reco_events = None
        self.reco_row = None
        self.create_widgets()
        self.load_reco()
        self.initialize_widget_values()
        ## Bubble detections are stored in the reco directory
        self.detected_bubbles = ReadBinary.ReadBlock(os.path.join(self.reco_directory, "ImageAnalysis_all.bin"))
        self.reset_event()

    def initialize_widget_values(self):
        values = sorted(self.reco_events.dtype.names) if self.reco_events is not None else ('')
        self.add_display_var_combobox['values'] = values
        self.manual_ped_config_directory.insert(0, self.config_file_directory)
        self.ped_config_file_path_combobox['values'] = self.get_configs()
        if os.path.isfile(self.ped_config_file_path_var):
            self.ped_config_file_path_combobox.insert(0, self.ped_config_file_path_var)
        if not os.path.exists(self.scan_directory):
            self.do_handscan_checkbutton['state'] = tk.DISABLED

    # reads config file and sets given values, otherwise sets default values
    def load_config_values(self, path):
        values = []
        for num in range(13):
            values.insert(num, None)
        defaults = []
        defaults.insert(0, str(self.raw_directory))
        defaults.insert(1, str(self.scan_directory))
        defaults.insert(2, str(self.reco_directory))
        defaults.insert(3, 'T3')
        defaults.insert(4, 'Images')
        defaults.insert(5, 0)
        defaults.insert(6, 4)
        defaults.insert(7, 3)
        defaults.insert(8, '30')
        defaults.insert(9, '50')
        defaults.insert(10, '70')
        defaults.insert(11, 'Piezo7')
        defaults.insert(12, 'Dytran')
        if os.path.isfile(path):
            set_default_vals = False
            f = open(path)
            counter = 0
            for line in f:
                parsed_line = line.rstrip('\n')
                if counter % 2 == 1:
                    values[int((counter - 1) / 2)] = parsed_line
                counter += 1
        else:
            logger.error('Config file not found at: ' + self.ped_config_file_path_var + '\nDefault values set')
        for num in range(len(values)):
            if values[num] is None:
                values[num] = defaults[num]
        self.raw_directory = values[0]
        self.scan_directory = values[1]
        self.reco_directory = values[2]
        self.plc_temp_var = values[3]
        self.images_relative_path = values[4]
        self.image_naming_convention = self.image_naming_conventions[int(values[5])]
        self.num_cams = int(values[6])
        self.image_orientation = self.image_orientations[int(values[7])]
        self.first_frame = values[8]
        self.init_frame = values[9]
        self.last_frame = values[10]
        self.piezo = values[11]
        self.dytran = values[12]
        self.frame = self.init_frame

    def reset_event(self):
        self.reco_row = None
        # Event of 20170719_0 Event 0
        self.row_index = 14177
        self.increment_event(1)

    def get_raw_events(self):
        try:
            # look for npy file in the custome directory
            # has to run convert_raw_to_npy.py first to generate the file

            self.raw_events = np.load(os.path.join(self.npy_directory,'raw_events.npy'))

        except FileNotFoundError:
            # this error should be handled when it crops up in the code
            raise FileNotFoundError

    def load_event_text(self):
        path = os.path.join(self.raw_directory, self.run, str(self.event), 'Event.txt')
        try:
            with open(path) as file:
                entries = file.readline().split()
                self.run_type = entries[2]
                self.run_type_label.set('run_type: ' + self.run_type)
                self.pset_label.set('pset: {:.1f}'.format(float(entries[9])))
                self.te_label.set('te: {:.1f}'.format(float(entries[10])))
        except FileNotFoundError:
            self.error += 'cannot find Event.txt\n'

    def load_plc_text(self):
        path = os.path.join(self.raw_directory, self.run, str(self.event), 'PLClog.txt')
        try:
            fields = linecache.getline(path, 2).split()
            index = fields.index(self.plc_temp_var)
            entries = linecache.getline(path, 7)
            self.temp_label.set(self.plc_temp_var + ': {:.1f}'.format(float(entries[index])))
        except ValueError:
            self.error += 'cannot find ' + self.plc_temp_var + ' in PLC log file\n'

    def reset_images(self):
        self.load_event_text()
        self.load_plc_text()
        self.load_fastDAQ_piezo(0)
        self.load_fastDAQ_piezo(1)
        self.load_PMT_traces()
        self.frame = self.init_frame
        self.diff_checkbutton_var.set(False)
        self.invert_checkbutton_var.set(False)
        for canvas in self.canvases:
            self.reset_zoom(canvas)
        #self.draw_fastDAQ_piezo_PMT_time()
        self.update_images()

    def get_image_path(self, cam, frame):
        # John Gresl 8/29/17: Set cam=0 since SBC has only one camera. Can comment out cam=0 in the event
        # that more cameras are added.
        cam=0
        if self.image_naming_convention == self.image_naming_conventions[0]:
            path = os.path.join(self.image_directory, 'cam{}_image{}.png'.format(cam, frame))
        elif self.image_naming_convention == self.image_naming_conventions[1]:
            # handle the leading spaces in the image names
            frame = '{:>3}'.format(frame)
            path = os.path.join(self.image_directory, 'cam{}image{}.bmp'.format(cam, frame))
        else:
            path = os.path.join(self.image_directory, 'cam{}_image{}.png'.format(cam, frame))
            self.error += ('Image naming convention not found\n')
        return path

    def update_images(self):
        error = ' '
        for canvas in self.canvases:
            path = self.get_image_path(canvas.cam, self.frame)
            image = self.load_image(path, canvas)
            zoom = '{:.1f}'.format(canvas.image_width / self.native_image_width)
            if self.diff_checkbutton_var.get():
                path = self.get_image_path(canvas.cam, self.first_frame)
                first_frame = self.load_image(path, canvas)
                image = ImageChops.difference(first_frame, image)
                template = 'frame: {} zoom: {}x (diff wrt {})                  {}/{}'
                bottom_text = template.format(self.frame, zoom, self.first_frame, self.run, self.event)
            else:
                template = 'frame: {} zoom: {}x                                   {}/{}'
                bottom_text = template.format(self.frame, zoom, self.run, self.event)
            if self.invert_checkbutton_var.get():
                image = ImageOps.invert(image)
            canvas.photo = ImageTk.PhotoImage(image)
            canvas.itemconfig(canvas.image, image=canvas.photo)
            canvas.itemconfig(canvas.bottom_text, text=bottom_text)
        self.draw_crosshairs()

    def load_image(self, path, canvas):
        try:
            image = Image.open(path)
            if self.image_orientation == '90':
                image = image.transpose(Image.ROTATE_90)
            if self.image_orientation == '180':
                image = image.transpose(Image.ROTATE_180)
            if self.image_orientation == '270':
                image = image.transpose(Image.ROTATE_270)
            self.trig_frame = 10
        except FileNotFoundError:
            self.trig_frame = -1
            self.error += ('did not find image at {}'.format(path) + '\n')
            image = Image.open(os.path.join(self.ped_directory, 'notfound.jpeg'))
        if self.show_bubbles_var.get():
            detected_bubbles = self.return_detected_bubbles()
            draw = ImageDraw.Draw(image)
            for bubble in detected_bubbles:
                draw.ellipse(list(bubble), outline=230)
        self.native_image_width, self.native_image_height = image.size
        image = image.resize((int(canvas.image_width), int(canvas.image_height)), self.antialias_checkbutton_var.get())
        image = image.crop((canvas.crop_left, canvas.crop_bottom, canvas.crop_right, canvas.crop_top))
        return image

    def load_run(self, run, event):
        if run == self.run and event == self.event:
            logger.info('no action taken (run and event are unchanged)')
        else:
            if len(np.argwhere((self.raw_events['run'] == run) & (self.raw_events['ev'] == event))) == 0:
                logger.error('invalid request: run {}, event {} does not exist'.format(run, event))
                self.update_run_entry()
                return
            logger.info('going to run {}, event {}'.format(run, event))
            self.run = run
            self.event = event
            if self.selected_events is None:
                self.row_index = self.get_row(self.raw_events)
            else:
                try:
                    self.row_index = self.get_row(self.selected_events)
                except IndexError:
                    logger.error('disabling cuts: requested run does not satisfy them')
                    self.selected_events = None
                    self.row_index = self.get_row(self.raw_events)
            self.update_run_entry()
            self.load_reco_row()
            self.image_directory = os.path.join(self.raw_directory, run, str(event), self.images_relative_path)
            self.reset_images()

    def add_display_var(self, var):
        if (self.reco_events is not None) and (var not in self.reco_events.dtype.names):
            logger.error('requested variable not in reco data: ' + var)
            return
        if var in [label['text'] for label, text, value in self.display_vars]:
            return
        label = tk.Label(self.bottom_frame_2, text=var)
        label.grid(row=len(self.display_vars) + 2, column=0)
        text = tk.StringVar(value=self.reco_row[var]) if self.reco_row else tk.StringVar(value='N/A')
        value = tk.Label(self.bottom_frame_2, textvariable=text, width=8)
        value.grid(row=len(self.display_vars) + 2, column=1, sticky='W')
        self.display_vars.append((label, text, value))

    def add_cut(self):
        field = ttk.Combobox(self.bottom_frame_1, width=3, values=sorted(self.reco_events.dtype.names))
        field.insert(0, 'nbub')
        field.grid(row=7 + len(self.cuts), column=0, columnspan=2, sticky='WE')
        operator = ttk.Combobox(self.bottom_frame_1, width=3, values=('>', '>=', '==', '<=', '<', '!='))
        operator.insert(0, '>=')
        operator.grid(row=7 + len(self.cuts), column=2, sticky='WE')
        value = tk.Entry(self.bottom_frame_1, width=5)
        value.insert(0, '0')
        value.grid(row=7 + len(self.cuts), column=3, sticky='WE')
        self.cuts.append((field, operator, value))

    def remove_cut(self):
        if not self.cuts:
            return
        for widget in self.cuts.pop():
            widget.destroy()
        self.apply_cuts()

    def remove_all_cuts(self):
        if not self.cuts:
            return
        while self.cuts:
            for widget in self.cuts.pop():
                widget.destroy()
        self.selected_events = None
        self.apply_cuts()

    def reset_cuts(self):
        for field, operator, value in self.cuts:
            field.delete(0, tk.END)
            operator.delete(0, tk.END)
            value.delete(0, tk.END)
        self.selected_events = None
        self.apply_cuts()

    def apply_cuts(self):
        if self.reco_events is None:
            logger.error('cannot apply cuts, reco data not found')
            return
        selection = []
        for field, operator, value in self.cuts:
            if field.get() == '' and operator.get() == '' and value.get() == '':
                continue
            if field.get() not in self.reco_events.dtype.names:
                logger.error('requested variable not in reco data')
                field.delete(0, tk.END)
                return
            dtype = self.reco_events[field.get()].dtype.str
            selection.append('(self.reco_events["{}"] {} {})'.format(
                field.get(),
                operator.get(),
                repr(value.get()) if 'U' in dtype else value.get()))  # add quotes if field datatype is string
        if len(selection) > 0:
            exec('self.selected_events = self.reco_events[{}]'.format(' & '.join(selection)))
            _, unique_rows = np.unique(self.selected_events[['run', 'ev']], return_index=True)
            self.selected_events = self.selected_events[unique_rows]  # get rid of multiple nbub entries
            if len(self.selected_events) == 0:
                logger.error('no events pass cuts')
                self.reset_cuts()
                return
            row = self.get_row(self.raw_events)
            try:
                events_left = self.raw_events[['run', 'ev']][row:]
                run, event = np.intersect1d(self.selected_events[['run', 'ev']], events_left)[0]
            except IndexError:
                logger.error('reached final event: starting over')
                run, event = self.selected_events[['run', 'ev']][0]
            self.run = run
            self.event = event
            self.reco_row = None
            self.row_index = self.get_row(self.selected_events) - 1
            self.increment_event(1)
        else:
            self.selected_events = None
            # self.row_index = self.get_row(self.raw_events)
            self.row_index = 0

    def get_row(self, events):
        return np.argwhere((events['run'] == self.run) & (events['ev'] == self.event)).ravel()[0]

    def increment_event(self, step):
        self.error = ''
        if self.selected_events is None:
            events = self.raw_events
        else:
            events = self.selected_events
        if (self.row_index + step) < 0:
            logger.error('reached first event: stopping here')
            self.reset_event()
            return
        if (self.row_index + step) >= len(events):
            logger.error('reached final event: starting over')
            self.reset_event()
            return
        self.row_index += step
        self.run = events[self.row_index]['run']
        self.event = events[self.row_index]['ev']
        self.update_run_entry()
        self.image_directory = os.path.join(self.raw_directory, self.run, str(self.event),
                                            self.images_relative_path)
        self.load_reco_row()
        self.n_PMT_trig.set(0)
        self.reset_images()
        self.goto_trig_frame()
        # printing out all error messages for an event at once
        event_path = "{raw}/{run}/{event}/Event.txt".format(raw=self.raw_directory,
                                                            run=self.run,
                                                            event=self.event)
        with open(event_path, "r") as event_txt:
            event_string = "Output from Event.txt:\n" + event_txt.read().strip()
        self.event_info_var.set(event_string)
        if not self.error == '':
            logger.error('This event had the following errors:\n' + self.error)

    ##find all dataset folders in the base directory of the form X*(size)-XX(year)-data
    def get_datasets(self):
        files = os.listdir(self.base_directory)
        i = 0
        while i < files.__len__():
            fileRegex = re.compile('\\w*-\\w*-data')
            if not fileRegex.match(files[i]):
                files.remove(files[i])
                i = i - 1
            i = i + 1
        return files

    # Returns a list of all config files in the config directory
    def get_configs(self):

        try: all_files = os.listdir(self.config_file_directory)
        except PermissionError: return None
        files = []
        for file in all_files:
            fileRegex = re.compile('\\w*-\\w*-ped_config.txt')
            if fileRegex.match(file):
                files.append(os.path.join(self.config_file_directory, str(file)))
        return files

    def update_run_entry(self):
        self.run_entry.delete(0, tk.END)
        self.run_entry.insert(0, self.run)
        self.event_entry.delete(0, tk.END)
        self.event_entry.insert(0, self.event)

    def set_init_dataset(self):
        counter = 0
        for dataset in self.get_datasets():
            if self.raw_directory.endswith(dataset):
                self.dataset_select.current(counter)
            counter += 1

    # Method for changing datasets from the Configuration tab
    def update_dataset(self):
        dataset = self.dataset_select.get()
        self.remove_all_cuts()
        if self.run_entry['state'] == tk.DISABLED:
            for child in self.bottom_frame_1.winfo_children():
                child.config(state=tk.NORMAL)
        try:
            split, end = dataset.split('-d', 1)
            values = self.ped_config_file_path_combobox['values']
            updated = False
            for value in values:
                if value.endswith(split + '-ped_config.txt'):
                    self.ped_config_file_path_combobox.delete(0, tk.END)
                    self.ped_config_file_path_combobox.insert(0, value)
                    updated = True
                    break
            if updated:
                self.config_combobox_update()
            ## if a config file is not found, this process will just set default values, error, and allow user to enter values
            else:
                self.raw_directory = os.path.join(self.base_directory, dataset)
                self.reco_directory = os.path.join(self.raw_directory, 'output')
                self.scan_directory = os.path.join(self.base_directory, 'scan_output_' + dataset[:-5])
                self.raw_directory_entry.delete(0, tk.END)
                self.raw_directory_entry.insert(0, self.raw_directory)
                self.reco_directory_entry.delete(0, tk.END)
                self.reco_directory_entry.insert(0, self.reco_directory)
                self.npy_directory_entry.delete(0, tk.END)
                self.npy_directory_entry.insert(0, self.npy_directory)
                self.scan_directory_entry.delete(0, tk.END)
                self.scan_directory_entry.insert(0, self.scan_directory)
                if os.path.exists(self.scan_directory):
                    self.do_handscan_checkbutton['state'] = tk.NORMAL
                else:
                    self.do_handscan_checkbutton['state'] = tk.DISABLED
                self.ped_config_file_path_combobox.delete(0, tk.END)
                self.ped_config_file_path_var = os.path.join(self.config_file_directory, split + '-ped_config.txt')
                if os.path.isfile(self.ped_config_file_path_var):
                    self.load_config_values(self.ped_config_file_path_var)
                else:
                    logger.error('Config file not found at: ' + self.ped_config_file_path_var +
                                 '\nDefault values set, please manually change these values if necessary')
                self.plc_temp_config_entry.delete(0, tk.END)
                self.plc_temp_config_entry.insert(0, self.plc_temp_var)
                self.relative_path_to_images_config_entry.delete(0, tk.END)
                self.relative_path_to_images_config_entry.insert(0, self.images_relative_path)
                self.image_naming_convention_select.delete(0, tk.END)
                self.image_naming_convention_select.insert(0, self.image_naming_convention)
                self.num_cams_config_entry.delete(0, tk.END)
                self.num_cams_config_entry.insert(0, self.num_cams)
                self.image_orientation_select.delete(0, tk.END)
                self.image_orientation_select.insert(0, self.image_orientation)
                self.first_frame_config_entry.delete(0, tk.END)
                self.first_frame_config_entry.insert(0, self.first_frame)
                self.init_frame_config_entry.delete(0, tk.END)
                self.init_frame_config_entry.insert(0, self.init_frame)
                self.last_frame_config_entry.delete(0, tk.END)
                self.last_frame_config_entry.insert(0, self.last_frame)
                self.update_num_cams()
                self.load_reco()
                self.get_raw_events()
                self.reset_event()
                self.set_init_dataset()
                self.reset_images()
        except FileNotFoundError:
            logger.error('Raw events not found for this dataset. Please ensure that the raw_events.npy file is present')
            self.num_cams = 0
            self.update_num_cams()
            if self.load_fastDAQ_piezo_checkbutton_vars[0].get():
                self.load_fastDAQ_piezo_checkbutton_vars[0].set(False)
                self.load_fastDAQ_piezo(0)
            if self.load_fastDAQ_piezo_checkbutton_vars[1].get():
                self.load_fastDAQ_piezo_checkbutton_vars[1].set(False)
                self.load_fastDAQ_piezo(1)
            for child in self.bottom_frame_1.winfo_children():
                child.config(state=tk.DISABLED)

    # Method for changing data directories
    def update_directories(self):
        try:
            self.raw_directory = self.raw_directory_entry.get()
            self.scan_directory = self.scan_directory_entry.get()
            self.reco_directory = self.reco_directory_entry.get()
            self.npy_directory = self.npy_directory_entry.get()
            if self.reco_directory_label['state'] != tk.DISABLED:
                if self.reco_version_combobox.get() == 'devel':
                    self.reco_directory = self.reco_directory.replace('current', 'devel')
                else:
                    self.reco_directory = self.reco_directory.replace('devel', 'current')
                self.reco_directory_entry.delete(0, tk.END)
                self.reco_directory_entry.insert(0, self.reco_directory)
            if not os.path.exists(self.scan_directory):
                self.do_handscan_checkbutton['state'] = tk.DISABLED
            else:
                self.do_handscan_checkbutton['state'] = tk.NORMAL
            if (not os.path.exists(self.raw_directory)) or (not os.path.exists(self.reco_directory)):
                raise FileNotFoundError
            self.get_raw_events()
            self.reset_event()
            self.remove_all_cuts()
            self.load_reco()
            self.set_init_dataset()
            self.num_cams = int(self.num_cams_config_entry.get())
            self.update_num_cams()
            self.reset_images()
            if self.run_entry['state'] == tk.DISABLED:
                for child in self.bottom_frame_1.winfo_children():
                    child.config(state=tk.NORMAL)
        except FileNotFoundError:
            logger.error(
                'One or more directories not found at given paths.\nPlease check paths, and also ensure that reco_events.npy and raw_events.npy '
                'exist in their respective directories')
            self.num_cams = 0
            self.update_num_cams()
            if self.load_fastDAQ_piezo_checkbutton_vars[0].get():
                self.load_fastDAQ_piezo_checkbutton_vars[0].set(False)
                self.load_fastDAQ_piezo(0)
            if self.load_fastDAQ_piezo_checkbutton_vars[1].get():
                self.load_fastDAQ_piezo_checkbutton_vars[1].set(False)
                self.load_fastDAQ_piezo(1)
            for child in self.bottom_frame_1.winfo_children():
                child.config(state=tk.DISABLED)

    # for when manual config path is updated
    def new_config_update(self):
        if os.path.exists(self.manual_ped_config_directory.get()):
            self.config_file_directory = self.manual_ped_config_directory.get()
            self.ped_config_file_path_combobox['values'] = self.get_configs()
        else:
            logger.error('Given config directory not found')

    # for when the config file path is changed
    def config_combobox_update(self):
        self.ped_config_file_path_var = self.ped_config_file_path_combobox.get()
        self.remove_all_cuts()
        self.load_config_values(self.ped_config_file_path_var)
        # update all of the widget values
        self.raw_directory_entry.delete(0, tk.END)
        self.raw_directory_entry.insert(0, self.raw_directory)
        self.scan_directory_entry.delete(0, tk.END)
        self.scan_directory_entry.insert(0, self.scan_directory)
        self.reco_directory_entry.delete(0, tk.END)
        self.reco_directory_entry.insert(0, self.reco_directory)
        self.npy_directory_entry.delete(0, tk.END)
        self.npy_directory_entry.insert(0, self.npy_directory)        
        if os.path.exists(self.scan_directory):
            self.do_handscan_checkbutton['state'] = tk.NORMAL
        else:
            self.do_handscan_checkbutton['state'] = tk.DISABLED
        self.plc_temp_config_entry.delete(0, tk.END)
        self.plc_temp_config_entry.insert(0, self.plc_temp_var)
        self.relative_path_to_images_config_entry.delete(0, tk.END)
        self.relative_path_to_images_config_entry.insert(0, self.images_relative_path)
        self.image_naming_convention_select.delete(0, tk.END)
        self.image_naming_convention_select.insert(0, self.image_naming_convention)
        self.num_cams_config_entry.delete(0, tk.END)
        self.num_cams_config_entry.insert(0, self.num_cams)
        self.image_orientation_select.delete(0, tk.END)
        self.image_orientation_select.insert(0, self.image_orientation)
        self.first_frame_config_entry.delete(0, tk.END)
        self.first_frame_config_entry.insert(0, self.first_frame)
        self.init_frame_config_entry.delete(0, tk.END)
        self.init_frame_config_entry.insert(0, self.init_frame)
        self.last_frame_config_entry.delete(0, tk.END)
        self.last_frame_config_entry.insert(0, self.last_frame)
        self.update_num_cams()
        self.load_reco()
        self.get_raw_events()
        self.reset_event()
        self.remove_all_cuts()
        self.set_init_dataset()
        self.reset_images()

    # for when values are changed wihtout updating config paths
    def update_vars_config(self):
        # plc var
        self.plc_temp_var = self.plc_temp_config_entry.get()
        self.load_plc_text()
        # image related vars
        self.image_orientation = self.image_orientation_select.get()
        self.num_cams = int(self.num_cams_config_entry.get())
        self.images_relative_path = self.relative_path_to_images_config_entry.get()
        self.image_naming_convention = self.image_naming_convention_select.get()
        self.init_frame = self.init_frame_config_entry.get()
        self.first_frame = self.first_frame_config_entry.get()
        self.last_frame = self.last_frame_config_entry.get()
        self.update_num_cams()
        self.reset_images()
        self.remove_all_cuts()
        self.reset_event()

    def update_num_cams(self):
        # reset the number of canvases
        for canvas in self.canvases:
            canvas.delete('all')
        self.canvases = []
        for cam in range(0, self.num_cams):
            canvas = tk.Canvas(self.camera_tab, width=self.init_image_width, height=self.init_image_height)
            canvas.bind('<ButtonPress-1>', self.on_button_press)
            canvas.zoom = 0
            canvas.image = canvas.create_image(0, 0, anchor=tk.NW, image=None)
            canvas.bottom_text = canvas.create_text(10, self.init_image_height - 25, anchor=tk.NW, text='', fill='red', font=DEFAULT_FONT)
            canvas.grid(row=0, column=1 * cam, columnspan=1, sticky='NW')
            canvas.cam = cam
            self.canvases.append(canvas)

    def load_reco_row(self, ibub=None):
        if self.reco_events is None:
            return
        if self.selected_events is not None:
            self.reco_row = self.selected_events[self.row_index]
        else:
            reco_index = self.raw_events[self.row_index]['reco index']
            if reco_index > 0:
                self.reco_row = self.reco_events[reco_index]
            else:
                self.reco_row = None
                self.toggle_reco_widgets(state=tk.DISABLED)
                for _, text, _ in self.display_vars:
                    text.set('N/A')
                return
        if ibub:
            offset = ibub - 1 if ibub > 1 else 0
            row = self.get_row(self.reco_events)
            self.reco_row = self.reco_events[row + offset]
        self.toggle_reco_widgets(state=tk.NORMAL)
        for label, text, _ in self.display_vars:
            var = label['text']
            dtype = self.reco_row[var].dtype.str
            text.set('{:.4f}'.format(self.reco_row[var]) if 'f' in dtype else self.reco_row[var])

    def toggle_reco_widgets(self, state):
        self.draw_crosshairs_button.config(state=state)
        for child in self.bottom_frame_2.winfo_children():
            child.config(state=state)

    def load_frame(self, frame):
        self.frame = str(frame)
        path = self.get_image_path(0, self.frame)
        if not os.path.isfile(path):
            self.frame = self.init_frame
        self.update_images()

    def goto_trig_frame(self):
        if self.trig_frame == -1:
            return
        date, run = self.run.split("_")
        date = np.int32(date)
        run = np.int32(run)
        event = np.int32(self.event)
        trig_frame = 10
        for n, (runid, ev, fr) in enumerate(
                zip(self.detected_bubbles["runid"], self.detected_bubbles["ev"], self.detected_bubbles["frame"])):
            if runid[0] == date and runid[1] == run and ev == event:
                trig_frame = self.detected_bubbles["trigFrame"][n]
        self.load_frame(trig_frame)
        return

    def reset_zoom(self, canvas):
        canvas.zoom = 0
        canvas.crop_left = 0
        canvas.crop_bottom = 0
        canvas.crop_right = self.init_image_width
        canvas.crop_top = self.init_image_height
        canvas.image_width = self.init_image_width
        canvas.image_height = self.init_image_height

    def on_button_press(self, event):
        canvas = event.widget
        old_width = canvas.image_width
        old_height = canvas.image_height
        canvas.zoom += 1
        if canvas.zoom > self.max_zoom:
            self.reset_zoom(canvas)
        else:
            canvas.image_width = 2 ** (canvas.zoom - 1) * self.native_image_width
            canvas.image_height = 2 ** (canvas.zoom - 1) * self.native_image_height
            if self.native_image_width < self.init_image_width:
                canvas.image_width = 2 ** (canvas.zoom) * self.init_image_width
            if self.native_image_height < self.init_image_height:
                canvas.image_height = 2 ** (canvas.zoom) * self.init_image_height
            new_center_x = (event.x + canvas.crop_left) * (canvas.image_width / old_width)
            new_center_y = (event.y + canvas.crop_bottom) * (canvas.image_height / old_height)
            if new_center_x < self.init_image_width / 2:
                # click was too far left, not enough new image for center to be here
                new_center_x = self.init_image_width / 2
            if new_center_x + self.init_image_width / 2 > canvas.image_width:
                # click was too far right, not enough new image for center to be here
                new_center_x = canvas.image_width - self.init_image_width / 2
            if new_center_y < self.init_image_height / 2:
                # click was too far up, not enough new image for center to be here
                new_center_y = self.init_image_height / 2
            if new_center_y + self.init_image_height / 2 > canvas.image_height:
                # click was too far down, not enough new image for center to be here
                new_center_y = canvas.image_height - self.init_image_height / 2
            canvas.crop_left = new_center_x - self.init_image_width / 2
            canvas.crop_bottom = new_center_y - self.init_image_height / 2
            canvas.crop_right = new_center_x + self.init_image_width / 2
            canvas.crop_top = new_center_y + self.init_image_height / 2
        self.update_images()

    def draw_crosshairs(self):
        for canvas in self.canvases:
            canvas.delete('crosshair')
        if not self.draw_crosshairs_var.get() or not self.reco_row:  # no reco row means we don't have reco data
            return
        if self.reco_row['nbub'] < 1:
            return
        for ibub in range(1, self.reco_row['nbub'] + 1):
            self.load_reco_row(ibub)
            for canvas in self.canvases:
                x_zoom = canvas.image_width / self.native_image_width
                y_zoom = canvas.image_height / self.native_image_height
                bubble_x = self.reco_row['hori{}'.format(canvas.cam)]
                bubble_y = self.reco_row['vert{}'.format(canvas.cam)]
                x = canvas.image_width - (bubble_x + canvas.crop_left / x_zoom) * x_zoom
                y = (bubble_y - canvas.crop_bottom / y_zoom) * y_zoom
                if self.image_orientation == '0':
                    x = (bubble_x - canvas.crop_left / x_zoom) * x_zoom
                    y = canvas.image_height - (bubble_y + canvas.crop_bottom / y_zoom) * y_zoom
                canvas.create_line(x - 11, y, x - 5, y, fill='red', tag='crosshair')
                canvas.create_line(x + 5, y, x + 11, y, fill='red', tag='crosshair')
                canvas.create_line(x, y - 11, x, y - 5, fill='red', tag='crosshair')
                canvas.create_line(x, y + 5, x, y + 11, fill='red', tag='crosshair')
                canvas.create_oval(x - 8, y - 8, x + 8, y + 8, outline='red', tag='crosshair')

    @staticmethod
    def toggle_widget(wid):
        # Toggles a tk widet from tk.NORMAL to tk.DISABLED or vice versa.
        wid["state"] = tk.NORMAL if wid["state"] == tk.DISABLED else tk.DISABLED
        return

    def toggle_widgets(self, fr):
        # Cycles through all widgets in a (fr)ame and toggles their state from tk.NORMAL to tk.DISABLED or vice versa.
        for child in fr.grid_slaves():
            self.toggle_widget(child)
        return

    def change_nbub(self):
        if self.nbub_button_var.get() > 1:
            for button in self.source_buttons:
                button.config(state=tk.ACTIVE)
            self.source_button_var.set(0)

    def load_fastDAQ_piezo(self, index):
        # Loads fastDAQ_piezo information. Variable index refers to either loading the top or bottom frames
        # All of the piezo variables should be stored in a list where the first element corresponds to the top,
        # and the second element corresponds to the bottom.
        if not self.load_fastDAQ_piezo_checkbutton_vars[index].get():
            self.destroy_children(self.piezo_tab_rights[index])
        path = os.path.join(self.raw_directory, self.run)
        self.fastDAQ_event = GetEvent(path, self.event, "fastDAQ")
        self.refresh_fastDAQ_piezo_choices()
        if self.load_fastDAQ_piezo_checkbutton_vars[index].get():
            self.toggle_widgets(self.piezo_checkbutton_frames[index])
        self.draw_fastDAQ_piezo(index)
        return

    def draw_all_fastDAQ_piezo_PMT_time(self):
        self.draw_fastDAQ_piezo_PMT_time(0)
        self.draw_fastDAQ_piezo_PMT_time(1)
        return

    def draw_fastDAQ_piezo_PMT_time(self, index):
        if self.piezo_ax[index] is None:
            return
        if not self.draw_time_in_fastDAQ_tab_var.get():
            try:
                self.piezo_line[index].remove()
                self.piezo_fig[index].canvas.draw()
            except:
                pass
            return
        if self.load_fastDAQ_piezo_checkbutton_vars[index].get():
            old_ylim = self.piezo_ax[index].get_ylim()
            if self.piezo_line[index] is not None:
                try:
                    self.piezo_line[index].remove()
                except:
                    pass
            self.piezo_line[index], = self.piezo_ax[index].plot([self.t0, self.t0], [-10, 10], "r-")
            self.piezo_ax[index].set_ylim(old_ylim)
            self.piezo_fig[index].canvas.draw()
        return

    def draw_fastDAQ_piezo(self, index):

        if not self.load_fastDAQ_piezo_checkbutton_vars[index].get():
            return
        if int(self.run_type) == 10:
            self.error += "Not allowed to view piezo data for run_type=10\n"
        self.piezos = self.get_active_piezo_checkboxes(index)
        try:
            self.piezo_cutoff_low = int(self.piezo_cutoff_low_entries[index].get())
            self.piezo_cutoff_high = int(self.piezo_cutoff_high_entries[index].get())
        except ValueError:
            logger.error(
                "Invalid types for cutoffs. Frequency cutoffs must be int, time cutoffs must be int or float.")
            return
        self.draw_filtered_piezo_trace(self.piezos, self.piezo_cutoff_low, self.piezo_cutoff_high, index)
        self.draw_fastDAQ_piezo_PMT_time(index)
        return

    def draw_filtered_piezo_trace(self, piezos, lowf, highf, index):
        board = 0
        try:
            if self.piezo_ax[index] is not None:
                for line in self.piezo_ax[index].lines[:]:
                    line.remove()
            for piezo in piezos:
                if piezo not in self.fastDAQ_event["fastDAQ"]["multiboards"][0] \
                        and piezo in self.fastDAQ_event["fastDAQ"]["multiboards"][1]:
                    board=1
                piezo_v = self.fastDAQ_event["fastDAQ"]["multiboards"][board][piezo]
                piezo_time = self.fastDAQ_event["fastDAQ"]["multiboards"][board]["time"]
                b, a = scipy.signal.butter(3, highf/len(piezo_v), "low")
                filtered_piezo_v = scipy.signal.lfilter(b, a, piezo_v)
                b, a = scipy.signal.butter(3, lowf/len(piezo_v), "high")
                filtered_piezo_v = scipy.signal.lfilter(b, a, filtered_piezo_v)
                if self.piezo_fig[index] is None:
                    self.piezo_fig[index], self.piezo_ax[index] = plt.subplots(figsize=(8, 3), dpi=100)
                self.piezo_ax[index].set_title(piezo)
                self.piezo_ax[index].set_xlabel("[s]")
                self.piezo_ax[index].set_ylabel("Amplitude [get units later]")
                self.piezo_ax[index].set_xlim(piezo_time[0], piezo_time[-1])
                plot_color="b"
                for color, cb in zip(self.piezo_colors, self.piezo_checkbuttons[index]):
                    if cb["text"] == piezo:
                        plot_color = color
                self.piezo_ax[index].plot(piezo_time, filtered_piezo_v, color=plot_color, label=piezo) # TODO: COLOR
                self.piezo_ax[index].set_ylim([-max(abs(filtered_piezo_v))*1.3, max(abs(filtered_piezo_v))*1.3])
                temp_legend=[]
                if self.piezo_plot_t0_checkbutton_vars[index].get():
                    #TODO: MAKE SURE THIS BUTTON IS DISABLED IF PMTs AREN'T LOADED
                    ### NEW VVV ###
                    # 1. Load binary reco file
                    acoustic_data_path = "{reco}/{run}/AcousticAnalysis_{run}.bin".format(reco=self.reco_directory,
                                                                                          run=self.run)
                    print("DEBUG:", acoustic_data_path)
                    print("DEBUG: FETCHING DATA")
                    self.acoustic_data = ReadBinary.ReadBlock(acoustic_data_path,max_file_size=2000)
                    print("DEBUG: DATA FETCHED")
                    self.piezo_t0 = self.acoustic_data["bubble_t0"][self.event]
                    print("DEBUG:", self.piezo_t0)
                    self.piezo_ax[index].axvline(x=self.piezo_t0[0], linestyle="dashed", color="r", label="t0[0]")
                    self.piezo_ax[index].axvline(x=self.piezo_t0[1], linestyle="dashed", color="b", label="t0[1]")
                    temp_legend = ["t0[0]", "t0[1]"]
                    ### NEW ^^^ ###

                    #TODO: Look at bubble t0 finding analysis to find out where it's going wrong
                    

                    ### OLD VVV ###
                    # if self.reco_row:
                    #     self.piezo_ax[index].axvline(x=self.reco_row["fastDAQ_t0"], linestyle="dashed", color="r", label="t0")
                    #     self.incremented_piezo_event = True
                    # else:
                    #     if self.incremented_piezo_event:
                    #         self.error += "t0 unavailable: No reco data found for current event. Piezo trace not drawn."
                    #     else:
                    #         logger.error("t0 unavailable: No reco data found for current event. Piezo trace not drawn.")
                    #     self.piezo_plot_t0_checkbutton_vars[index].set(False)
                    #     self.incremented_piezo_event = False
                    ### OLD ^^^ ###



                temp_sticky = tk.NW if index == 0 else tk.SW
                self.place_graph_and_toolbar(figure=self.piezo_fig[index], master=self.piezo_tab_rights[index],
                                             sticky=temp_sticky)
            # MAKE THE LEGEND
            if self.piezo_ax[index] is not None:
                self.piezo_ax[index].legend(bbox_to_anchor=(0.8, 0.3), loc='upper center')

        except (KeyError, IndexError):

            self.error += "Piezo data not found.\n"
            logger.error("Piezo data not found.")
            self.destroy_children(self.piezo_tab_rights[index])
            canvas = tk.Canvas(self.piezo_tab_rights[index], width=self.init_image_width, height=self.init_image_height)
            self.load_fastDAQ_piezo_checkbuttons[index].toggle()
            self.reset_zoom(canvas)
        return

    def get_active_piezo_checkboxes(self, index):
        out = [ ]
        n_cb = len(self.piezo_checkbox_vars[index])
        for n in range(n_cb):
            if self.piezo_checkbox_vars[index][n].get():
                out.append(self.piezo_checkbuttons[index][n]["text"])
        #out.append(self.piezo_checkbuttons[index][n]["text"] if self.piezo_checkbox_vars[index][n].get() else 0 for n in range(len(self.piezo_checkbox_vars[index])))
        return out

    def draw_piezos_from_checkbuttons(self, index=0):
        # Go through the list of checkbuttons for a given index and plot them on the index's graph
        # with the correct color from self.piezo_colors .. or something.
        n_cb = len(self.piezo_checkbox_vars[index])
        for n in range(n_cb):
            if self.piezo_checkbox_vars[index][n].get():
                label = self.piezo_checkbuttons[index][n]["text"]
                print("(Not) Plotting", label)


        return


    def refresh_fastDAQ_piezo_choices(self):
        if self.fastDAQ_event["fastDAQ"]["loaded"] is False:
            return
        board = 0
        new_choices = list(self.fastDAQ_event["fastDAQ"]["multiboards"][board].keys())
        exclude = ["time", "loaded", "bindata", "caldata", "multiboards"]
        choices = [choice for choice in new_choices if choice not in exclude]
        del self.piezo_checkbox_vars[0][:]
        del self.piezo_checkbox_vars[1][:]
        del self.piezo_checkbuttons[0][:]
        del self.piezo_checkbuttons[1][:]
        for n in range(len(choices)):
            self.piezo_checkbox_vars[0].append(tk.BooleanVar(master=self.piezo_checkbutton_frames[0], value=0))
            self.piezo_checkbuttons[0].append(tk.Checkbutton(master=self.piezo_checkbutton_frames[0],
                                                             text=choices[n], variable=self.piezo_checkbox_vars[0][-1],
                                                             command=lambda : self.draw_fastDAQ_piezo(0),
                                                             state = tk.DISABLED))
            self.piezo_checkbuttons[0][-1].grid(row=n, column=0, sticky=tk.N)

            self.piezo_checkbox_vars[1].append(tk.BooleanVar(master=self.piezo_checkbutton_frames[1], value=0))
            self.piezo_checkbuttons[1].append(tk.Checkbutton(master=self.piezo_checkbutton_frames[1],
                                                             text=choices[n], variable=self.piezo_checkbox_vars[1][-1],
                                                             command=lambda : self.draw_fastDAQ_piezo(1),
                                                             state=tk.DISABLED))
            self.piezo_checkbuttons[1][-1].grid(row=n, column=0, sticky=tk.N)
        return

    def draw_time_in_piezo(self):
        # Checks if a line has already been draw, deletes it, and draws another. This is to save from having to
        # redraw the piezo graphs every time we move forward a PMT trigger.
        return

    def increment_PMT_trigger(self, n):
        self.n_PMT_trig.set(value=self.n_PMT_trig.get() + n)
        self.load_PMT_traces()
        self.draw_fastDAQ_piezo_PMT_time(0)
        self.draw_fastDAQ_piezo_PMT_time(1)
        return

    def forward_PMT_trigger(self):
        # TODO: what happens when you reach the final trigger?
        # TODO: what happens when you reach the final trigger?
        self.increment_PMT_trigger(1)
        return

    def back_PMT_trigger(self):
        self.increment_PMT_trigger(-1)
        return


    def jump_to_t0_trigger(self):
        # Scan through all the times and jump to the PMT trigger closest to t=0
        if not self.draw_pmt_traces_var.get():
            return
        if self.PMT_data is not None and self.align_data is not None:
            n_max = self.PMT_data["t0_sec"][:, 0].shape[0]
            times = [ ]
            align_t0_sec = self.align_data["PMT_trigt0_sec"][self.event]
            align_t0_frac = self.align_data["PMT_trigt0_frac"][self.event]
            for n in range(n_max):
                # Improve this by using a binary search method in the future.
                # Do this now as a proof of concept
                trace_t0_sec = self.PMT_data["t0_sec"][n, 0]
                trace_t0_frac = self.PMT_data["t0_frac"][n, 0]
                times.append(abs(self.return_pmt_time((trace_t0_sec, trace_t0_frac), (align_t0_sec, align_t0_frac))))
            min_index = times.index(min(times))
            self.increment_PMT_trigger(min_index-self.n_PMT_trig.get())
        return


    @staticmethod
    def disable_widgets(widget_list):
        for widget in widget_list:
            widget["state"] = tk.DISABLED
        return

    @staticmethod
    def enable_widgets(widget_list):
        for widget in widget_list:
            widget["state"] = tk.NORMAL
        return

    def load_PMT_traces(self):
        self.destroy_children(self.pmt_graph_frame)
        if not self.draw_pmt_traces_var.get():
            self.draw_time_in_fastDAQ_tab_var.set(0)
            self.piezo_plot_t0_checkbutton_var_top.set(0)
            self.piezo_plot_t0_checkbutton_var_bottom.set(0)
            self.disable_widgets(self.bottom_frame_3_bottom.grid_slaves() + [self.piezo_plot_t0_checkbutton_top,self.piezo_plot_t0_checkbutton_bottom])
            return
        self.enable_widgets(self.bottom_frame_3_bottom.grid_slaves() + [self.piezo_plot_t0_checkbutton_top,
                                                                        self.piezo_plot_t0_checkbutton_bottom])
        pmt_data_path = "{raw}/{run}/{event}/PMTtraces.bin".format(raw=self.raw_directory,
                                                                   run=self.run, event=self.event)
        align_data_path = "{reco}/{run}/PMTfastDAQalignment_{run}.bin".format(reco=self.reco_directory,
                                                                              run=self.run)
        if self.current_PMT_event != self.event:
            try:
                self.PMT_data = ReadBinary.ReadBlock(pmt_data_path, max_file_size=2000) # Create some kind of flag to check if the file is
                                                                    # currently loaded into memory so we can
                self.jump_to_t0_trigger_button["state"] = tk.NORMAL
            except:
                raise Exception("Cannot get PMT data") #TODO: Raise custom exception to be caught by draw_PMT_traces
            try:
                self.align_data = ReadBinary.ReadBlock(align_data_path, max_file_size=2000)
            except:
                raise Exception("Cannot find align data") # TODO: Same ^
        if self.n_PMT_trig.get()< 0:
            logger.error("PMT trigger cannot be less than 0.")
            self.n_PMT_trig.set(0)
            return
        if self.n_PMT_trig.get() > self.PMT_data["t0_sec"][:,0].shape[0]:
            logger.error("PMT trigger cannot be greater than", str(self.PMT_data["t0_sec"][:,0].shape[2]))
            self.n_PMT_trig.set(self.PMT_data["t0_sec"][:,0].shape[0])
            return
        self.draw_PMT_traces(self.PMT_data, self.align_data)
        self.current_PMT_event = self.event
        return

    def draw_PMT_traces(self, PMT_data, align_data):
        n = self.n_PMT_trig.get()
        f = matplotlib.pyplot.Figure(figsize=(11, 2.5), dpi=100)
        a = f.add_subplot(111)
        if PMT_data is None:
            logger.error("PMT_data not found.")
            return
        if align_data is None:
            logger.error("PMT alignment data not found.")
        xd = np.arange(PMT_data["traces"].shape[2]) * PMT_data["dt"][n, 0]
        yd_fine = PMT_data["traces"][n, 0, :] * PMT_data["v_scale"][n, 0] + \
                  PMT_data["v_offset"][n, 0]
        # yd_coarse = PMT_data["traces"][n, 1, :] * PMT_data["v_scale"][n, 1] + \
        #             PMT_data["v_offset"][n, 1]
        if align_data is not None:
            trace_t0_sec = PMT_data["t0_sec"][n, 0]
            trace_t0_frac = PMT_data["t0_frac"][n, 0]
            align_t0_sec = align_data["PMT_trigt0_sec"][self.event]
            align_t0_frac = align_data["PMT_trigt0_frac"][self.event]
            # Sanity Check
            if align_t0_frac < 0:
                logger.error("PMT fastDAQ alignment failed. Unable to recover calibrated time of PMT trigger.")
                self.t0 = -1.0
            else:
                self.t0 = self.return_pmt_time((trace_t0_sec, trace_t0_frac), (align_t0_sec, align_t0_frac))
        a.plot(xd, yd_fine)
        a.set_xlabel("[s]")
        a.set_ylabel("PMT Amplitude [mV]")
        a.set_title("PMT Trace at t0 = {}".format(self.t0))
        self.place_graph_and_toolbar(figure=f, master=self.pmt_graph_frame)
        return

    @staticmethod
    def return_pmt_time(trace_t0, align_t0):
        # Returns the calibrated time of the PMT trace. Both variables are 2-tuples, with the first element being the
        # "seconds" component, and the 2nd element being the "fractional" component.
        return (trace_t0[0] - align_t0[0]) + (trace_t0[1] - align_t0[1])


    @staticmethod
    def place_graph_and_toolbar(figure, master, sticky=tk.NW):
        # Given a figure and a tkinter master where the graph will go, this method
        # places and grids the graph and a useful toolbar for saving, zooming, etc.
        canvas = FigureCanvasTkAgg(figure, master=master)
        tb_frame = tk.Frame(master=master)
        tb_frame.grid(row=1, column=0, sticky=sticky)
        tb = CustomToolbar(canvas, tb_frame)
        tb.update()
        canvas._tkcanvas.grid(row=0, column=0, sticky=tk.N)
        return

    def destroy_children(self, frame):
        for widget in frame.winfo_children():
            widget.destroy()
        return

    @staticmethod
    def center_point_to_corner_points(x,y,rad):
        # Given a circle's center coordinate (x,y) and radius, this method returns
        # the upper left coordinate (x1, y1) and the lower right coordinate (x2,y2)
        # in the form of a 4x1 numpy float array [x1, y1, x2, y2].
        return np.array([x-rad, y-rad, x+rad, y+rad], dtype=float)

    def return_detected_bubbles(self):
        # John Gresl 8/30/2017. (Eventually) Returns the locations of any bubbles detected by an
        # external routine in the form of an array of bubble locations. For example, if two bubbles
        # are detected, this method will return [ [xa1, ya1, xa2, ya2], [xb1, yb1, xa2, yb2] ]
        # where (x1,y1) is the top left coordinate, (x2,y2) is the bottom right coordinate, and a
        # and b represent the first and second bubble, respectively.
        bubbles = [ ]
        date, run = self.run.split("_")
        date = np.int32(date)
        run = np.int32(run)
        event = np.int32(self.event)
        frame = np.int32(self.frame)
        for n, (runid, ev, fr) in enumerate(zip(self.detected_bubbles["runid"], self.detected_bubbles["ev"], self.detected_bubbles["frame"])):
            if runid[0] == date and runid[1] == run and ev == event and fr == frame:
                ## Add bubble locations to the list
                x_coord = self.detected_bubbles["ipix"][n]
                y_coord = self.detected_bubbles["jpix"][n]
                bubbles.append(self.center_point_to_corner_points(x_coord, y_coord, 18))
        return bubbles

    def load_reco(self):
        self.reco_row = None
        self.reco_events = None
        path = os.path.join(self.npy_directory, 'reco_eventsNOTYET.npy')
        if not os.path.isfile(path):
            logger.error('cannot find reco_data.npy, reco data will be disabled')
            self.toggle_reco_widgets(state=tk.DISABLED)
            for _, text, _ in self.display_vars:
                text.set('N/A')
            return
        logger.info('using reco data from {}'.format(path))
        events = np.load(path)
        if len(events) == 0:
            logger.error('could not find raw data for any reco events')
            return
        self.reco_events = events

    def do_handscan(self):
        if self.do_handscan_checkbutton_var.get():
            file = 'scan_{}_{}_{}.txt'.format(self.run, getpass.getuser(), time.strftime('%a_%b_%d_%H_%M_%S_%Y'))
            self.scan_file = os.path.join(self.scan_directory, file)
            state = tk.NORMAL
            self.scanner_name.set('scanner: ' + getpass.getuser())
        else:
            if os.path.isfile(self.scan_file) and (os.stat(self.scan_file).st_size == 0):
                os.remove(self.scan_file)
            state = tk.DISABLED
            self.scanner_name.set('')
        for i in range(0, 8):
            self.nbub_button[i].config(state=state)
        for button in self.source_buttons:
            button.config(state=state)
        self.isgoodtrigger_button.config(state=state)
        self.crosshairsgood_button.config(state=state)
        self.comment_label.config(state=state)
        self.comment_entry.config(state=state)
        self.submit_scan_button.config(state=state)

    def submit_scan(self):
        if ((str(self.scanner_name.get())[9:] == '') or
                (str(self.source_button_var.get()) == '-1') or
                (str(self.nbub_button_var.get()) == '-1')):
            logger.error('did not complete scan selections')
            return

        with(open(self.scan_file, 'a+')) as file:
            file.seek(0)
            if not file.read(1):
                file.write('Output from ped hand scanning.\n')
                file.write(
                    'run  ev  scanner  scan_source  scan_nbub  scan_trigger  scan_crosshairsgood  scan_comment\n')
                file.write('%s  %d  %s  %d  %d  %d  %d  %s\n1\n\n\n')
            file.write(self.run + '  ' +
                       str(self.event) + '  ' +
                       str(self.scanner_name.get())[9:] + '  ' +
                       str(self.source_button_var.get()) + ' ' +
                       str(self.nbub_button_var.get()) + ' ' +
                       str(int(self.isgoodtrigger_checkbutton_var.get())) + ' ' +
                       str(int(self.crosshairsgood_checkbutton_var.get())) + ' ' +
                       '\'' + self.comment_entry.get() + '\'\n')
        self.nbub_button_var.set(-1)
        self.source_button_var.set(-1)
        self.isgoodtrigger_checkbutton_var.set(True)
        self.crosshairsgood_checkbutton_var.set(True)
        self.comment_entry.delete(0, tk.END)
        self.comment_entry.insert(0, '')
        self.increment_event(1)

    def create_widgets(self):
        self.notebook = ttk.Notebook(self, padding=[0, 0, 0, 0])
        self.camera_tab = tk.Frame(self.notebook)
        self.piezo_tab = tk.Frame(self.notebook)
        self.PMT_tab = tk.Frame(self.notebook)
        self.config_tab = tk.Frame(self.notebook)

        self.notebook.add(self.camera_tab, text='camera')
        self.notebook.add(self.piezo_tab, text='piezo')
        self.notebook.add(self.PMT_tab, text="PMT")
        self.notebook.add(self.config_tab, text='Configuration')
        self.notebook.grid(row=0, column=0, columnspan=5)

        # Setup frames to be used on the top (in tabs)
        # Cameras tab
        self.canvases = []
        for cam in range(0, self.num_cams):
            canvas = tk.Canvas(self.camera_tab, width=self.init_image_width, height=self.init_image_height)
            canvas.bind('<ButtonPress-1>', self.on_button_press)
            canvas.zoom = 0

            canvas.image = canvas.create_image(0, 0, anchor=tk.NW, image=None)
            canvas.bottom_text = canvas.create_text(10, self.init_image_height - 25, anchor=tk.NW, text='', fill='red', font=DEFAULT_FONT)
            canvas.grid(row=0, column=1 * cam, columnspan=1, sticky='NW')
            canvas.cam = cam
            self.canvases.append(canvas)
        # Piezos tab
        # First setup frames for piezos tab
        self.piezo_tab_left = tk.Frame(self.piezo_tab, bd=5, relief=tk.SUNKEN)
        self.piezo_tab_left.grid(row=0, column=0, sticky='NW')

        self.piezo_tab_left_top = tk.Frame(self.piezo_tab_left, bd=2, relief=tk.SUNKEN)
        self.piezo_tab_left_top.grid(row=0, column=0, sticky=tk.N)
        self.piezo_tab_left_bottom = tk.Frame(self.piezo_tab_left, bd=2, relief=tk.SUNKEN)
        self.piezo_tab_left_bottom.grid(row=1, column=0, sticky=tk.N, pady=(140, 0))

        self.piezo_tab_right = tk.Frame(self.piezo_tab, bd=5, relief=tk.SUNKEN)
        self.piezo_tab_right.grid(row=0, column=1, sticky='NW')

        self.piezo_tab_right_top = tk.Frame(self.piezo_tab_right, bd=2, relief=tk.SUNKEN)
        self.piezo_tab_right_top.grid(row=0, column=0, sticky=tk.N)
        self.piezo_tab_right_bottom = tk.Frame(self.piezo_tab_right, bd=2, relief=tk.SUNKEN)
        self.piezo_tab_right_bottom.grid(row=1, column=0, sticky=tk.S)
        self.piezo_tab_rights = [self.piezo_tab_right_top,
                                 self.piezo_tab_right_bottom]
        self.piezo_ax = [None, None]
        self.piezo_line = [None, None]
        self.piezo_fig = [None, None]
        # Now within the piezos frames setup stuff
        self.load_fastDAQ_piezo_checkbutton_top = tk.Checkbutton(
            self.piezo_tab_left_top,
            text='Load fastDAQ',
            variable=self.load_fastDAQ_piezo_checkbutton_var_top, #top
            command=lambda: self.load_fastDAQ_piezo(0))
        self.load_fastDAQ_piezo_checkbutton_top.grid(row=0, column=0, columnspan=2, sticky='WE')
        self.load_fastDAQ_piezo_checkbutton_bottom = tk.Checkbutton(
            self.piezo_tab_left_bottom,
            text="Load fastDAQ",
            variable=self.load_fastDAQ_piezo_checkbutton_var_bottom,
            command=lambda: self.load_fastDAQ_piezo(1))
        self.load_fastDAQ_piezo_checkbutton_bottom.grid(row=0, column=0, columnspan=2, sticky=tk.EW)
        self.load_fastDAQ_piezo_checkbuttons = [self.load_fastDAQ_piezo_checkbutton_top,
                                                self.load_fastDAQ_piezo_checkbutton_bottom]
        self.piezo_label_top = tk.Label(self.piezo_tab_left_top, text='Piezo:')
        self.piezo_label_top.grid(row=1, column=0, sticky='WE')

        self.piezo_label_bottom = tk.Label(self.piezo_tab_left_bottom, text='Piezo:')
        self.piezo_label_bottom.grid(row=1, column=0, sticky='WE')
        self.piezo_labels = [self.piezo_label_top,
                             self.piezo_label_bottom]


        self.piezo_checkbutton_frame_top = tk.Frame(master=self.piezo_tab_left_top)
        self.piezo_checkbutton_frame_top.grid(row=1, column=1, sticky=tk.N)
        self.piezo_checkbutton_frame_bot = tk.Frame(master=self.piezo_tab_left_bottom)
        self.piezo_checkbutton_frame_bot.grid(row=1, column=1, sticky=tk.N)
        self.piezo_checkbutton_frames = [self.piezo_checkbutton_frame_top,
                                         self.piezo_checkbutton_frame_bot]
        choices = ["Dytran"]
        self.piezo_checkbox_vars_top = [ ]
        self.piezo_checkbox_vars_bot = [ ]
        self.piezo_checkbuttons_top = [ ]
        self.piezo_checkbuttons_bot = [ ]
        self.piezo_colors = ['b','g','c','k','y','m']*2
        self.piezo_checkbox_vars = [self.piezo_checkbox_vars_top,
                                    self.piezo_checkbox_vars_bot]
        self.piezo_checkbuttons = [self.piezo_checkbuttons_top,
                                   self.piezo_checkbuttons_bot]
        # Build the checkbuttons



        #self.piezo_dropdown_var_top = tk.StringVar(master=self.piezo_tab_right_top, value = choices[0])
        self.piezo_dropdown_var_bottom = tk.StringVar(master=self.piezo_tab_right_top, value=choices[0])
        self.piezo_dropdown_vars = [None,
                                    self.piezo_dropdown_var_bottom]

        #self.piezo_dropdown_top = tk.OptionMenu(master=self.piezo_tab_left_top,
        #                                    value=choices,
        #                                    variable=self.piezo_dropdown_var_top)
        #self.piezo_dropdown_top.grid(row=1, column=1, sticky=tk.EW)
        #self.piezo_dropdown_bottom = tk.OptionMenu(master=self.piezo_tab_left_bottom,
        #                                        value=choices,
        #                                        variable=self.piezo_dropdown_var_bottom)
        #self.piezo_dropdown_bottom.grid(row=1, column=1, sticky=tk.EW)
        #self.piezo_dropdowns = [None,
        #                        self.piezo_dropdown_bottom]

        self.piezo_cutoff_low_label_top = tk.Label(self.piezo_tab_left_top, text='Freq cutoff low:')
        self.piezo_cutoff_low_label_top.grid(row=2, column=0, sticky='WE')
        self.piezo_cutoff_low_label_bottom = tk.Label(self.piezo_tab_left_bottom, text='Freq cutoff low:')
        self.piezo_cutoff_low_label_bottom.grid(row=2, column=0, sticky='WE')
        self.piezo_cutoff_low_labels = [self.piezo_cutoff_low_label_top,
                                        self.piezo_cutoff_low_label_bottom]

        self.piezo_cutoff_low_entry_top = tk.Entry(self.piezo_tab_left_top, width=12)
        self.piezo_cutoff_low_entry_top.insert(0, self.piezo_cutoff_low)
        self.piezo_cutoff_low_entry_top.grid(row=2, column=1, sticky='WE')
        self.piezo_cutoff_low_entry_bottom = tk.Entry(self.piezo_tab_left_bottom, width=12)
        self.piezo_cutoff_low_entry_bottom.insert(0, self.piezo_cutoff_low)
        self.piezo_cutoff_low_entry_bottom.grid(row=2, column=1, sticky='WE')
        self.piezo_cutoff_low_entries = [self.piezo_cutoff_low_entry_top,
                                         self.piezo_cutoff_low_entry_bottom]

        self.piezo_cutoff_high_label_top = tk.Label(self.piezo_tab_left_top, text='Freq cutoff high:')
        self.piezo_cutoff_high_label_top.grid(row=3, column=0, sticky='WE')
        self.piezo_cutoff_high_label_bottom = tk.Label(self.piezo_tab_left_bottom, text='Freq cutoff high:')
        self.piezo_cutoff_high_label_bottom.grid(row=3, column=0, sticky='WE')
        self.piezo_cutoff_high_labels = [self.piezo_cutoff_high_label_top,
                                         self.piezo_cutoff_high_label_bottom]

        self.piezo_cutoff_high_entry_top = tk.Entry(self.piezo_tab_left_top, width=12)
        self.piezo_cutoff_high_entry_top.insert(0, self.piezo_cutoff_high)
        self.piezo_cutoff_high_entry_top.grid(row=3, column=1, sticky='WE')
        self.piezo_cutoff_high_entry_bottom = tk.Entry(self.piezo_tab_left_bottom, width=12)
        self.piezo_cutoff_high_entry_bottom.insert(0, self.piezo_cutoff_high)
        self.piezo_cutoff_high_entry_bottom.grid(row=3, column=1, sticky='WE')
        self.piezo_cutoff_high_entries = [self.piezo_cutoff_high_entry_top,
                                          self.piezo_cutoff_high_entry_bottom]
        # Commented out JG 10/12/2017. No need to have these buttons and options when we can just
        # chose our time window using the zoom tool. May want to revert, just uncomment this if you do.
        # self.piezo_timerange_checkbutton_top = tk.Checkbutton(
        #     self.piezo_tab_left_top, text='Full time window',
        #     variable=self.piezo_timerange_checkbutton_var_top,
        #     command=lambda: self.draw_fastDAQ_piezo(0))
        # self.piezo_timerange_checkbutton_top.grid(row=6, column=0, columnspan=2, sticky='WE')
        # self.piezo_timerange_checkbutton_bottom = tk.Checkbutton(
        #     self.piezo_tab_left_bottom, text='Full time window',
        #     variable=self.piezo_timerange_checkbutton_var_bottom,
        #     command=lambda: self.draw_fastDAQ_piezo(1))
        # self.piezo_timerange_checkbutton_bottom.grid(row=6, column=0, columnspan=2, sticky='WE')
        # self.piezo_timerange_checkbuttons = [self.piezo_timerange_checkbutton_top,
        #                                      self.piezo_timerange_checkbutton_bottom]
        #
        # self.piezo_beginning_time_label_top = tk.Label(self.piezo_tab_left_top, text='Beginning Time:')
        # self.piezo_beginning_time_label_top.grid(row=4, column=0, sticky='WE')
        # self.piezo_beginning_time_label_bottom = tk.Label(self.piezo_tab_left_bottom, text='Beginning Time:')
        # self.piezo_beginning_time_label_bottom.grid(row=4, column=0, sticky='WE')
        # self.piezo_beginning_time_labels = [self.piezo_beginning_time_label_top,
        #                                     self.piezo_beginning_time_label_bottom]
        #
        # self.piezo_beginning_time_entry_top = tk.Entry(self.piezo_tab_left_top, width=12)
        # self.piezo_beginning_time_entry_top.insert(0, self.piezo_beginning_time)
        # self.piezo_beginning_time_entry_top.grid(row=4, column=1, sticky='WE')
        # self.piezo_beginning_time_entry_bottom = tk.Entry(self.piezo_tab_left_bottom, width=12)
        # self.piezo_beginning_time_entry_bottom.insert(0, self.piezo_beginning_time)
        # self.piezo_beginning_time_entry_bottom.grid(row=4, column=1, sticky='WE')
        # self.piezo_beginning_time_entries = [self.piezo_beginning_time_entry_top,
        #                                      self.piezo_beginning_time_entry_bottom]
        #
        # self.piezo_ending_time_label_top = tk.Label(self.piezo_tab_left_top, text='Ending Time:')
        # self.piezo_ending_time_label_top.grid(row=5, column=0, sticky='WE')
        # self.piezo_ending_time_label_bottom = tk.Label(self.piezo_tab_left_bottom, text='Ending Time:')
        # self.piezo_ending_time_label_bottom.grid(row=5, column=0, sticky='WE')
        # self.piezo_ending_time_labels = [self.piezo_ending_time_label_top,
        #                                  self.piezo_ending_time_label_bottom]
        #
        # self.piezo_ending_time_entry_top = tk.Entry(self.piezo_tab_left_top, width=12)
        # self.piezo_ending_time_entry_top.insert(0, self.piezo_ending_time)
        # self.piezo_ending_time_entry_top.grid(row=5, column=1, sticky='WE')
        # self.piezo_ending_time_entry_bottom = tk.Entry(self.piezo_tab_left_bottom, width=12)
        # self.piezo_ending_time_entry_bottom.insert(0, self.piezo_ending_time)
        # self.piezo_ending_time_entry_bottom.grid(row=5, column=1, sticky='WE')
        # self.piezo_ending_time_entries = [self.piezo_ending_time_entry_top,
        #                                   self.piezo_ending_time_entry_bottom]

        self.piezo_plot_t0_checkbutton_top = tk.Checkbutton(
            self.piezo_tab_left_top,
            text='Show t0',
            variable=self.piezo_plot_t0_checkbutton_var_top,
            command=lambda: self.draw_fastDAQ_piezo(0))
        self.piezo_plot_t0_checkbutton_top.grid(row=7, column=0, columnspan=2, sticky='WE')
        self.piezo_plot_t0_checkbutton_bottom = tk.Checkbutton(
            self.piezo_tab_left_bottom,
            text='Show t0',
            variable=self.piezo_plot_t0_checkbutton_var_bottom,
            command=lambda: self.draw_fastDAQ_piezo(1))
        self.piezo_plot_t0_checkbutton_bottom.grid(row=7, column=0, columnspan=2, sticky='WE')
        self.piezo_plot_t0_checkbuttons = [self.piezo_plot_t0_checkbutton_top,
                                           self.piezo_plot_t0_checkbutton_bottom]

        self.reload_fastDAQ_piezo_button_top = tk.Button(self.piezo_tab_left_top, text='reload',
                                                     command=lambda: self.draw_fastDAQ_piezo(0))
        self.reload_fastDAQ_piezo_button_top.grid(row=8, column=0, sticky='WE')
        self.reload_fastDAQ_piezo_button_bottom = tk.Button(self.piezo_tab_left_bottom, text='reload',
                                                     command=lambda: self.draw_fastDAQ_piezo(1))
        self.reload_fastDAQ_piezo_button_bottom.grid(row=8, column=0, sticky='WE')
        self.reload_fastDAQ_piezo_buttons = [self.reload_fastDAQ_piezo_button_top,
                                             self.reload_fastDAQ_piezo_button_bottom]
        # PMT tab
        self.pmt_settings_frame = tk.Frame(master=self.PMT_tab, bd=5, relief=tk.SUNKEN)
        self.pmt_settings_frame.grid(row=0, column=0, sticky=tk.NW)
        self.pmt_graph_frame = tk.Frame(master=self.PMT_tab, bd=5, relief=tk.SUNKEN)
        self.pmt_graph_frame.grid(row=0, column=1, rowspan=2, sticky=tk.N)
        self.draw_pmt_traces_var = tk.BooleanVar(value=0)
        self.draw_pmt_traces_checkbutton = tk.Checkbutton(master=self.pmt_settings_frame,
                                                          text="Draw PMT Traces",
                                                          variable = self.draw_pmt_traces_var,
                                                          command = self.load_PMT_traces)
        self.draw_pmt_traces_checkbutton.grid(row=0, column=0, sticky=tk.N)

        # configuration tab
        # setup frames within tab
        self.config_tab_main = tk.Frame(self.config_tab, bd=5, relief=tk.SUNKEN)
        self.config_tab_main.grid(row=0, column=0, sticky='NW')

        self.config_tab_vars = tk.Frame(self.config_tab, bd=5, relief=tk.SUNKEN)
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

        # Setup frames to be used on the bottom
        self.bottom_frame_1 = tk.Frame(self, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_1.grid(row=1, column=0, sticky='NW')
        self.bottom_frame_2 = tk.Frame(self, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_2.grid(row=1, column=1, sticky='NW')
        self.bottom_frame_3 = tk.Frame(self, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_3.grid(row=1, column=2, sticky='NW')
        self.bottom_frame_3_top = tk.Frame(self.bottom_frame_3, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_3_top.grid(row=0, column=0, sticky=tk.N)
        self.bottom_frame_3_bottom = tk.Frame(self.bottom_frame_3, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_3_bottom.grid(row=1, column=0, sticky=tk.N)
        self.bottom_frame_4 = tk.Frame(self, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_4.grid(row=1, column=3, sticky='NW')
        self.bottom_frame_5 = tk.Frame(master=self.bottom_frame_1, bd=5, relief=tk.SUNKEN)
        self.bottom_frame_5.grid(row=10, column=0, columnspan=4, sticky=tk.N)

        self.run_label = tk.Label(self.bottom_frame_1, text='run:')
        self.run_label.grid(row=1, column=0, sticky='WE')

        self.run_entry = tk.Entry(self.bottom_frame_1, width=12)
        self.run_entry.grid(row=1, column=1, sticky='WE')

        self.event_label = tk.Label(self.bottom_frame_1, text='event:')
        self.event_label.grid(row=1, column=2, sticky='WE')

        self.event_entry = tk.Entry(self.bottom_frame_1, width=5)
        self.event_entry.grid(row=1, column=3, sticky='WE')

        self.go_button = tk.Button(self.bottom_frame_1, text='Go', command=self.load_run)
        self.go_button['command'] = lambda: self.load_run(self.run_entry.get(), int(self.event_entry.get()))
        self.go_button.grid(row=1, column=4, sticky='WE')

        self.back_event = tk.Button(self.bottom_frame_1, text='back event', command=lambda: self.increment_event(-1))
        self.back_event.grid(row=2, column=0, columnspan=2, sticky='WE')

        self.forward_event = tk.Button(self.bottom_frame_1, text='forward event',
                                       command=lambda: self.increment_event(1))
        self.forward_event.grid(row=2, column=2, columnspan=2, sticky='WE')

        self.back_1000events_button = tk.Button(self.bottom_frame_1, text='back 1000 events')
        self.back_1000events_button['command'] = lambda: self.increment_event(-1000)
        self.back_1000events_button.grid(row=3, column=0, columnspan=2, sticky='WE')

        self.forward_1000events_button = tk.Button(self.bottom_frame_1, text='forward 1000 events')
        self.forward_1000events_button['command'] = lambda: self.increment_event(1000)
        self.forward_1000events_button.grid(row=3, column=2, columnspan=2, sticky='WE')

        self.fill_run_type = tk.Label(self.bottom_frame_1, textvariable=self.run_type_label, width=11)
        self.fill_run_type.grid(row=4, column=0, sticky='WE')

        self.fill_pset = tk.Label(self.bottom_frame_1, textvariable=self.pset_label, width=10)
        self.fill_pset.grid(row=4, column=1, sticky='WE')

        self.fill_temp = tk.Label(self.bottom_frame_1, textvariable=self.temp_label, width=10)
        self.fill_temp.grid(row=4, column=2, sticky='WE')

        self.fill_te = tk.Label(self.bottom_frame_1, textvariable=self.te_label, width=10)
        self.fill_te.grid(row=4, column=3, sticky='WE')

        self.reset_cuts_button = tk.Button(self.bottom_frame_1, text='reset cuts', command=self.reset_cuts)
        self.reset_cuts_button.grid(row=6, column=0, columnspan=2, sticky='WE')

        self.apply_cuts_button = tk.Button(self.bottom_frame_1, text='apply cuts', command=self.apply_cuts)
        self.apply_cuts_button.grid(row=6, column=2, columnspan=2, sticky='WE')

        self.add_cut_button = tk.Button(self.bottom_frame_1, text='add cut', command=self.add_cut)
        self.add_cut_button.grid(row=5, column=0, columnspan=2, sticky='WE')

        self.remove_cut_button = tk.Button(self.bottom_frame_1, text='delete cut', command=self.remove_cut)
        self.remove_cut_button.grid(row=5, column=2, columnspan=2, sticky='WE')

        self.display_reco_label = tk.Label(self.bottom_frame_2, text='Variables from merged_all')
        self.display_reco_label.grid(row=0, column=0, sticky='WE')

        self.add_display_var_combobox = ttk.Combobox(self.bottom_frame_2)
        self.add_display_var_combobox.grid(row=1, column=0, sticky='WE')

        self.add_display_var_button = tk.Button(
            self.bottom_frame_2,
            text='add',
            command=lambda: self.add_display_var(self.add_display_var_combobox.get()))
        self.add_display_var_button.grid(row=1, column=1, sticky='WE')

        self.display_vars = []
        self.add_display_var('nbub')
        #         self.add_display_var('getBub_success')
        self.add_display_var('fastDAQ_t0')
        self.add_display_var('te')

        self.back_frame_button = tk.Button(self.bottom_frame_3_top, text='back frame')
        self.back_frame_button['command'] = lambda: self.load_frame(int(self.frame) - 1)
        self.back_frame_button.grid(row=0, column=0, sticky='WE')

        self.forward_frame_button = tk.Button(self.bottom_frame_3_top, text='forward frame')
        self.forward_frame_button['command'] = lambda: self.load_frame(int(self.frame) + 1)
        self.forward_frame_button.grid(row=0, column=1, sticky='WE')

        self.reset_images_button = tk.Button(self.bottom_frame_3_top, text='reset image', command=self.reset_images)
        self.reset_images_button.grid(row=0, column=2, sticky='WE')

        self.first_frame_button = tk.Button(self.bottom_frame_3_top, text='first frame')
        self.first_frame_button['command'] = lambda: self.load_frame(self.first_frame)
        self.first_frame_button.grid(row=1, column=0, sticky='WE')

        self.last_frame_button = tk.Button(self.bottom_frame_3_top, text='last frame')
        self.last_frame_button['command'] = lambda: self.load_frame(self.last_frame)
        self.last_frame_button.grid(row=1, column=1, sticky='WE')

        self.trig_frame_button = tk.Button(self.bottom_frame_3_top, text='trig frame')
        self.trig_frame_button["command"] = self.goto_trig_frame
        self.trig_frame_button.grid(row=1, column=2, sticky='WE')

        self.antialias_checkbutton = tk.Checkbutton(
            self.bottom_frame_3_top,
            text='antialias',
            variable=self.antialias_checkbutton_var,
            command=self.update_images)
        self.antialias_checkbutton.grid(row=2, column=0, sticky='WE')

        self.diff_checkbutton = tk.Checkbutton(
            self.bottom_frame_3_top,
            text='diff frame',
            variable=self.diff_checkbutton_var,
            command=self.update_images)
        self.diff_checkbutton.grid(row=2, column=1, sticky='WE')

        self.invert_checkbutton = tk.Checkbutton(
            self.bottom_frame_3_top,
            text='invert',
            variable=self.invert_checkbutton_var,
            command=self.update_images)
        self.invert_checkbutton.grid(row=2, column=2, sticky='WE')

        self.show_bubbles_checkbutton = tk.Checkbutton(
            master=self.bottom_frame_3_top,
            text="show bubbles",
            variable=self.show_bubbles_var,
            command=self.update_images
        )
        self.show_bubbles_checkbutton.grid(row=3, column=0, sticky=tk.W)

        self.draw_crosshairs_button = tk.Checkbutton(
            self.bottom_frame_3_top,
            text='draw crosshairs',
            variable=self.draw_crosshairs_var,
            command=self.draw_crosshairs,
            state=tk.DISABLED)
        self.draw_crosshairs_button.grid(row=3, column=1, sticky='E')
        self.pmt_trigger_label = tk.Label(master=self.bottom_frame_3_bottom,
                                          text="PMT Trigger")
        self.pmt_trigger_label.grid(row=0, column=0, sticky=tk.EW)
        self.pmt_trigger_entry = tk.Entry(master=self.bottom_frame_3_bottom,
                                          textvariable=self.n_PMT_trig)
        self.pmt_trigger_entry.grid(row=0, column=1, sticky=tk.EW)
        self.pmt_go_button = tk.Button(master=self.bottom_frame_3_bottom,
                                       text="Go",
                                       command=self.reset_images)
        self.pmt_go_button.grid(row=0, column=2, sticky=tk.N)
        self.draw_time_in_fastDAQ_tab_var = tk.BooleanVar(value=0)
        self.draw_time_in_fastDAQ_tab_checkbutton = tk.Checkbutton(master=self.bottom_frame_3_bottom,
                                                                   text="Draw PMT time in fastDAQ tab",
                                                                   variable=self.draw_time_in_fastDAQ_tab_var,
                                                                   command=self.draw_all_fastDAQ_piezo_PMT_time)
        self.draw_time_in_fastDAQ_tab_checkbutton.grid(row=1, column=0, columnspan=2, sticky=tk.W)
        self.back_PMT_trigger_button = tk.Button(master=self.bottom_frame_3_bottom,
                                                 text="Back PMT trigger",
                                                 command=self.back_PMT_trigger)
        self.back_PMT_trigger_button.grid(row=2, column=0, sticky=tk.EW)
        self.forward_PMT_trigger_button = tk.Button(master=self.bottom_frame_3_bottom,
                                                    text="Forward PMT trigger",
                                                    command=self.forward_PMT_trigger)
        self.forward_PMT_trigger_button.grid(row=2, column=1, sticky=tk.N)

        self.jump_to_t0_trigger_button = tk.Button(master=self.bottom_frame_3_bottom,
                                                   text="Jump to t~0 trigger.",
                                                   command=self.jump_to_t0_trigger,
                                                   state=tk.DISABLED)
        self.jump_to_t0_trigger_button.grid(row=3, column=0, sticky=tk.N)
        self.do_handscan_checkbutton = tk.Checkbutton(
            self.bottom_frame_4,
            text='do handscan',
            variable=self.do_handscan_checkbutton_var,
            command=self.do_handscan)
        self.do_handscan_checkbutton.grid(row=0, column=0, columnspan=2, sticky='WE')

        self.scanner_name_label = tk.Label(self.bottom_frame_4, textvariable=self.scanner_name)
        self.scanner_name_label.grid(row=0, column=2, columnspan=2, sticky='WE')

        self.nbub_button = []
        max_bub = 8
        for i, text in zip(range(0, max_bub), [' bubbles'] * (max_bub - 1) + ['+ bubbles']):
            self.nbub_button.append(
                tk.Radiobutton(
                    self.bottom_frame_4,
                    text=str(i) + text,
                    state=tk.DISABLED,
                    variable=self.nbub_button_var,
                    command=self.change_nbub,
                    value=i))
            self.nbub_button[i].grid(row=1 + i, column=0, columnspan=2, sticky='WE')

        self.source_buttons = []
        for i, text in enumerate(['bulk event', 'wall with tail', 'surface', 'collar', 'other']):
            button = tk.Radiobutton(self.bottom_frame_4, state=tk.DISABLED, variable=self.source_button_var, value=i)
            button['text'] = text
            button.grid(row=i + 1, column=2, columnspan=2, sticky='WE')
            self.source_buttons.append(button)

        self.isgoodtrigger_button = tk.Checkbutton(self.bottom_frame_4, variable=self.isgoodtrigger_checkbutton_var)
        self.isgoodtrigger_button['text'] = 'Is good trigger?'
        self.isgoodtrigger_button['state'] = tk.DISABLED
        self.isgoodtrigger_button.grid(row=7, column=2, columnspan=2, sticky='WE')

        self.crosshairsgood_button = tk.Checkbutton(self.bottom_frame_4, variable=self.crosshairsgood_checkbutton_var)
        self.crosshairsgood_button['text'] = 'Crosshairs good?'
        self.crosshairsgood_button['state'] = tk.DISABLED
        self.crosshairsgood_button.grid(row=8, column=2, columnspan=2, sticky='WE')

        self.comment_label = tk.Label(self.bottom_frame_4, text='Comment:', state=tk.DISABLED)
        self.comment_label.grid(row=9, column=0, sticky='WE')

        self.comment_entry = tk.Entry(self.bottom_frame_4, width=15, state=tk.DISABLED)
        self.comment_entry.insert(0, '')
        self.comment_entry.grid(row=9, column=1, columnspan=2, sticky='WE')

        self.submit_scan_button = tk.Button(self.bottom_frame_4, state=tk.DISABLED, text='Submit and Go')
        self.submit_scan_button['command'] = self.submit_scan
        self.submit_scan_button.grid(row=9, column=3, sticky='WE')

        self.event_info_var = tk.StringVar(master=self.bottom_frame_5, value="")
        self.event_info_label = tk.Label(master=self.bottom_frame_5, textvariable=self.event_info_var, width=40, wraplength=250)
        self.event_info_label.grid(row=0, column=0, sticky=tk.N)
ROOT = tk.Tk()
ROOT.lift()
ROOT.attributes('-topmost', True)
ROOT.after_idle(ROOT.attributes, '-topmost', False)
WIDTH = 1650  # width for the Tk root
HEIGHT = 950  # height for the Tk root
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
