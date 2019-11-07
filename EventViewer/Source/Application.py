import os
import tkinter as tk


class Application(tk.Frame):
    def __init__(self, master=None,):
        tk.Frame.__init__(self, master)
        self.grid()
        self.max_zoom = 3
        self.init_image_width = 400
        self.init_image_height = 625
        # TODO: What are these for? Maybe define in configs.py instead. Maybe zoom level too?
        self.native_image_width = None
        self.native_image_height = None
        self.error = "" # Errors are appended to this string.
        self.plc_temp_var = "T3"
        # TODO: Do we need this^?
        self.images_relative_path = "Images" # Relative to each event
        self.image_naming_convention = "cam0_image0.png"
        self.image_orientation = "0"
        self.num_cams = 4
        # TODO: Should this be 1?
        self.first_frame = "30"
        self.init_frame = "50"
        self.last_frame = "70"
        # TODO: What are these 3^?
        self.piezo = "Piezo7"
        # TODO: What is this^?
        self.ped_config_file_path_var = os.path.join(self.config_file_directory, "SBC-17-ped_config.txt")
        return