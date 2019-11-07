
import os
import tkinter as tk
from tkinter import filedialog
# from PIL import ImageTk, Image
from EventViewer.Source.WidgetTooltip import WidgetToolTip
from EventViewer.Source.ErrorWindow import ErrorWindow
from EventViewer.Source.GUIDefaults import *
import EventViewer.Source.Setup.nOptionDialogBox as nOptionDialogBox


class EvVarSetup(object):
    def __init__(self):
        self.ok_clicked = False
        self.overwriteQ = False
        self.root = tk.Tk()
        self.root.eval("tk::PlaceWindow %s center" %
                       self.root.winfo_pathname(self.root.winfo_id()))
        self.root.title("EventViewer Setup")
        self.setup_h = 1920
        self.setup_w = 1080
        self.canvas = tk.Canvas(master=self.root,
                                height=self.setup_h,
                                width=self.setup_w)
        self.canvas.grid(row=0,
                         column=0,
                         sticky=tk.NSEW)
        self.root.resizable(False, False)
        self.path_label = tk.Label(master=self.canvas,
                                   text="Paths",
                                   font=HEADER_FONT)
        self.path_label.grid(row=0,
                             column=0,
                             sticky=tk.E)

        self.help_button = tk.Label(master=self.canvas,
                                    text="?",
                                    font=HEADER_FONT,
                                    relief=tk.SUNKEN,
                                    width=2)
        WidgetToolTip(self.help_button,
                      text="Hover over an item for help. If you're unsure of one of the parameters " +
                           "it's best to leave it blank.")
        self.help_button.grid(row=0,
                              column=2,
                              sticky=tk.N)

        self.path_labels = [" Raw Directory", "Scan Directory", "Reco Directory"]
        self.path_tooltips = ["Directory where raw_events.npy is located",
                              "Directory where handscans will go. Disabled currently",
                              "Directory of reconstructed data (where reco_events.npy is located)"]
        self.disabled_paths = [False, True, False]  # Corresponding to the order in path_pathels

        self.path_stringvars = []
        for n, label in enumerate(self.path_labels):
            t_frame = tk.Frame(master=self.canvas)
            t_frame.grid(row=n + 1,
                         column=1,
                         sticky=tk.NSEW)
            lab = tk.Label(master=t_frame,
                           text=label + ":",
                           font=DEFAULT_FONT)
            lab.grid(row=0,
                     column=0,
                     sticky=tk.W)
            WidgetToolTip(lab,
                          self.path_tooltips[n])
            self.path_stringvars.append(tk.StringVar(master=t_frame))
            tk.Entry(master=t_frame,
                     width=50,
                     textvariable=self.path_stringvars[n],
                     font=DEFAULT_FONT).grid(row=0,
                                             column=1,
                                             sticky=tk.NSEW)
            tk.Button(master=t_frame,
                      text="Browse...",
                      font=DEFAULT_FONT,
                      command=lambda i=n: self.on_browse(self.path_stringvars[i]))\
                            .grid(row=0,
                                  column=2,
                                  sticky=tk.E)

        # DISABLE SCAN DIRECTORY FOR NOW (Note: .grid_slaves() returns in chronologically reversed grid order
        for n in range(len(self.disabled_paths)):
            rev_n = len(self.disabled_paths)-n-1  # grid_slaves is in reversed order
            if self.disabled_paths[n]:
                for grid_slave in self.canvas.grid_slaves()[rev_n].grid_slaves():
                    grid_slave.config(state=tk.DISABLED)
        # DISABLE SCAN DIRECTORY FOR NOW (Note: .grid_slaves() returns in chronologically reversed grid order

        self.parameter_label = tk.Label(master=self.canvas,
                                        text="Parameters",
                                        font=HEADER_FONT)
        self.parameter_label.grid(row=len(self.path_labels) + 2,
                                  column=0,
                                  sticky=tk.E)
        self.parameter_labels = ["PLC temperature var", "Relative Path to images", "Image naming convention",
                                 "Number of cameras", "Image Orientation", "First frame", "Trigger frame",
                                 "Piezo", "Dytran"]
        self.parameter_defaults = ["T3", " ", "1", "2", "0", "1", "10", "Piezo7", "Dytran"]
        self.parameter_tooltips = ["...",
                                   "Leave this blank for now",
                                   "0 for 'cam0_image0.png', 1 for 'cam0image  0.bmp'",
                                   "Number of cameras",
                                   "0 for '0 deg', 1 for '90 deg', 2 for '180 deg', 3 for '270 deg'",
                                   "The initial frame you want EventViewer to display",
                                   "The trigger frame (defined by your DAQ)",
                                   "Piezo variable name",
                                   "Dytran variable name"]
        self.disabled_parameters = [False, True, False, False, False,
                                    False, False, False, False]
        self.parameter_stringvars = []
        self._max_width = max([len(p) + 1 for p in self.parameter_labels])
        for n, parameter in enumerate(self.parameter_labels):
            p_frame = tk.Frame(master=self.canvas)
            p_frame.grid(row=len(self.path_labels) + n + 3,
                         column=1,
                         sticky=tk.E)
            lab = tk.Label(master=p_frame,
                           text=" " * (self._max_width - len(parameter) - 1) + parameter + ":",
                           width=self._max_width, font=DEFAULT_FONT)
            WidgetToolTip(lab,
                          self.parameter_tooltips[n])
            lab.grid(row=0,
                     column=0,
                     sticky=tk.W)
            self.parameter_stringvars.append(tk.StringVar(master=p_frame,
                                                          value=self.parameter_defaults[n]))
            tk.Entry(master=p_frame,
                     width=50,
                     textvariable=self.parameter_stringvars[n],
                     font=DEFAULT_FONT).grid(row=0,
                                             column=1,
                                             sticky=tk.W)
        # DISABLE RELATIVE IMAGE PATH (Note: .grid_slaves() returns in chronologically reversed grid order
        for n in range(len(self.disabled_parameters)):
            rev_n = len(self.disabled_parameters)-n-1  # grid_slaves is in reversed order
            if self.disabled_parameters[n]:
                for grid_slave in self.canvas.grid_slaves()[rev_n].grid_slaves():
                    grid_slave.config(state=tk.DISABLED)
        # DISABLE RELATIVE IMAGE PATH (Note: .grid_slaves() returns in chronologically reversed grid order


        self.output_frame = tk.Frame(master=self.root)
        self.output_frame.grid(row=1,
                               column=0, sticky=tk.NSEW)
        self.output_label = tk.Label(master=self.output_frame, text=" Write config file to:", font=DEFAULT_FONT)
        self.output_label.grid(row=0, column=0, sticky=tk.NSEW)
        self.output_path_var = tk.StringVar(master=self.output_frame, value="/")
        self.output_path_entry = tk.Entry(master=self.output_frame,
                                          textvariable=self.output_path_var,
                                          font=DEFAULT_FONT,
                                          width=54)
        self.output_path_entry.grid(row=0, column=1, sticky=tk.NSEW)
        tk.Button(master=self.output_frame,
                  text="Browse...",
                  font=DEFAULT_FONT,
                  command=lambda : self.on_browse(self.output_path_var, askfor="file"))\
            .grid(row=0,
                  column=2,
                  sticky=tk.NSEW, padx=5)
        WidgetToolTip(self.output_label, "Place where Event-Viewer config file will go. "+
                                         "For default config, don't change this.")

        self.button_frame = tk.Frame(master=self.root)
        self.button_frame.grid(row=2,
                               column=0,
                               sticky=tk.E)
        self.cancel_button = tk.Button(master=self.button_frame,
                                       text="Cancel",
                                       font=DEFAULT_FONT,
                                       command=self.on_cancel)
        self.cancel_button.grid(row=0,
                                column=0,
                                sticky=tk.W)
        self.ok_button = tk.Button(master=self.button_frame,
                                   text="OK",
                                   font=DEFAULT_FONT,
                                   command=self.on_ok)
        self.ok_button.grid(row=0,
                            column=1,
                            sticky=tk.E)

        self.root.mainloop()
        self.file_paths = [self.path_stringvars[n].get() for n in range(len(self.path_labels))]


        return

    def overwrite(self, click_ev):
        self.overwriteQ = True

    def dont_overwrite(self, click_ev):
        self.overwriteQ = False

    def generate_config_str(self):
        template = '''Raw Directory:
{raw_dir}
Scan Directory:
{scan_dir}
Reco Directory:
{reco_dir}
PLC temperature var:
{plc_temp}
Relative Path to images:
{image_path}
Image naming convention (0 for 'cam0_image0.png', 1 for 'cam0image  0.bmp'):
{image_conv}
Number of Cameras:
{n_cam}
Image Orientation (0 for '0', 1 for '90', 2 for '180', 3 for '270'): 
{image_orient}
first frame: 
{first_frame}
trig frame: 
{trig_frame}
last frame: 
{last_frame}
Piezo: 
{piezo}
Dytran:
{dytran}
'''
        output = template.format(raw_dir=)

        return

    def on_ok(self):
        error = ""
        for n in range(len(self.path_labels)):
            if not self.disabled_paths[n]:  # ... Meaning the path is enabled.
                if not os.path.isdir(self.path_stringvars[n].get()):
                    error += "-{} = '{}' is not a valid directory.\n\n".\
                        format(self.path_labels[n].strip(),
                               self.path_stringvars[n].get())
        if os.path.isfile(self.output_path_var.get()):
            nOptionDialogBox.nOptionDialogBox(prompt="A file already exists at\n'{}'\nDo you wish to overwrite?"\
                                              .format(self.output_path_var.get()),
                                              options=["Yes", "No"],
                                              commands=[self.overwrite, self.dont_overwrite],
                                              master=self.root)

        if error != "":
            ErrorWindow(self.root, error)
            return


        self.ok_clicked = True

        # self.root.destroy()
        return

    def on_cancel(self):
        self.root.destroy()
        return

    def on_browse(self, textvar, askfor="dir"):
        if askfor == "dir":
            d = filedialog.askdirectory(initialdir="/")
        else:
            d = filedialog.asksaveasfilename(initialdir="/")
        textvar.set(d)
        return


if __name__ == "__main__":
    e = EvVarSetup()
