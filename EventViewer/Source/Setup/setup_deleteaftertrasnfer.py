## John Gresl
## j.gresl12@gmail.com


## The start of every new window should be centered in the screen. To do this we can use a bit
## of c-code which is already implemented in Tkinter and deals with multiply monitors correctly.
## So when instantiating a new .Tk() object, do this
# root = tk.Tk()
# root.eval("tk::PlaceWindow %s center" % root.winfo_pathname(root.winfo_id()))

import os
import tkinter as tk
from tkinter import filedialog
from PIL import ImageTk, Image
from EventViewer.Source.WidgetTooltip import WidgetToolTip
from EventViewer.Source.ErrorWindow import ErrorWindow



DEFAULT_FONT = ("Consolas", 14)
HEADER_FONT = (DEFAULT_FONT[0], int(DEFAULT_FONT[1]*1.3))
ROOT_TITLE = "EventViewer Setup"

def run_setup():
    exit_code, _ = run_intro()
    if exit_code != 1:
        return exit_code
    ask_and_convert_npy()
    idk_yet()


    return exit_code, 0

def check_paths(path_arr):
    return all([os.path.isdir(p) for p in path_arr])


def ask_and_convert_npy():
    global exit_code
    exit_code = 0
    def _on_yes():
        print("Yass")
        return
    def _on_no():
        global exit_code
        exit_code = 1
        root.destroy()
        return


    root = tk.Tk()
    root.eval("tk::PlaceWindow %s center" % root.winfo_pathname(root.winfo_id()))

    root.title(ROOT_TITLE)
    root.resizable(False, False)
    npy_h = 200; npy_w = 400
    canvas = tk.Canvas(master=root, height=npy_h, width=npy_w)
    _img = ImageTk.PhotoImage(Image.open("../Resources/placeholder.jpg"))
    canvas.create_image(0, 0, anchor=tk.NW, image=_img)
    canvas.pack()
    _text = "Would you like to create your raw\n"+\
            "and reco .npy files now? You only\n"+\
            "have to do this once per data set\n"+\
            "or change of your merged_all.bin.\n"

    textx, texty = rel_to_abs(0.5, 0.3, npy_w, npy_h)
    canvas.create_text(textx, texty, text=_text, anchor=tk.CENTER, font=DEFAULT_FONT)

    button_frame = tk.Frame(master=canvas)
    button_frame.place(relx=0.5, rely=0.7, anchor="center")

    no_button = tk.Button(master=button_frame, text="No", font=DEFAULT_FONT, width=5, command=_on_no)
    no_button.grid(row=0, column=0, sticky=tk.E)
    yes_button = tk.Button(master=button_frame, text="Yes", font=DEFAULT_FONT, width=5, command=_on_yes)
    yes_button.grid(row=0, column=1, sticky=tk.W)

    root.mainloop()
    return


def idk_yet():
    # We want to gather the data paths
    global exit_code
    exit_code = 0
    def _on_browse(n):
        dir = filedialog.askdirectory(initialdir="/")
        path_stringvars[n].set(dir)
        return
    def _on_ok():
        global exit_code
        if check_paths([path_stringvars[i].get() for i in [0,2]]): ## Skip scan dir (1) since that isn't setup yet
            print("Paths look ok.")
        else:
            ErrorWindow(root, "Paths are not ok.")
        return
    def _on_cancel():
        print("cancel")
        return
    root = tk.Tk()
    root.eval("tk::PlaceWindow %s center" % root.winfo_pathname(root.winfo_id()))
    root.title(ROOT_TITLE)
    idk_h = 1920; idk_w = 1080
    canvas = tk.Canvas(master=root, height=idk_h, width=idk_w)
    canvas.grid(row=0, column=0, sticky=tk.NSEW)
    root.resizable(False, False)
    path_label = tk.Label(master=canvas, text="Paths", font=HEADER_FONT)
    path_label.grid(row=0, column=0, sticky=tk.E)
    #_img = tk.PhotoImage(file="../Resources/help.gif")
    #_imgx, _imgy = rel_to_abs(0, 0, idk_w, idk_h)
    #canvas.create_image(_imgx, _imgy, image=_img)
    help_button = tk.Label(master=canvas, text="?", font=HEADER_FONT, relief=tk.SUNKEN, width=2)
    WidgetToolTip(help_button, "Hover over an item for help. If you're unsure of one of the parameters"+\
                     " it's best to leave it blank.")
    help_button.grid(row=0, column=2, sticky=tk.N)

    path_labels = [" Raw Directory", "Scan Directory", "Reco Directory"]
    path_tooltips = ["Directory where raw_events.npy is located",
                     "Directory where handscans will go?? Disabled currently",
                     "Directory of reconstructed data (where reco_events.npy is located)"]

    path_stringvars = []
    for n, label in enumerate(path_labels):
        t_frame = tk.Frame(master=canvas)
        t_frame.grid(row=n+1, column=1, sticky=tk.NSEW)
        lab = tk.Label(master=t_frame, text=label+":", font=DEFAULT_FONT)
        lab.grid(row=0, column=0, sticky=tk.W)
        WidgetToolTip(lab, path_tooltips[n])
        path_stringvars.append(tk.StringVar(master=t_frame))
        tk.Entry(master=t_frame, width=50, textvariable=path_stringvars[n], font=DEFAULT_FONT)\
            .grid(row=0, column=1, sticky=tk.NSEW)
        tk.Button(master=t_frame, text="Browse...", font=DEFAULT_FONT, command= lambda n=n: _on_browse(n))\
            .grid(row=0, column=2, sticky=tk.E)

    ## DISABLE SCAN DIRECTORY FOR NOW (Note: .grid_slaves() returns in chronologically reversed grid order
    canvas.grid_slaves()[1].grid_slaves()[0].config(state=tk.DISABLED)
    canvas.grid_slaves()[1].grid_slaves()[1].config(state=tk.DISABLED)
    canvas.grid_slaves()[1].grid_slaves()[2].config(state=tk.DISABLED)

    ## DISABLE SCAN DIRECTORY FOR NOW (Note: .grid_slaves() returns in chronologically reversed grid order

    parameter_label = tk.Label(master=canvas, text="Parameters", font=HEADER_FONT)
    parameter_label.grid(row=len(path_labels)+2, column=0, sticky=tk.E)
    parameter_labels = ["PLC temperature var", "Relative Path to images", "Image naming convention",
                        "Number of cameras", "Image Orientation", "First frame", "Trigger frame",
                        "Piezo", "Dytran"]
    parameter_defaults = ["T3", " ", "1", "2", "0", "1", "10", "Piezo7", "Dytran"]
    parameter_tooltips = ["...",
                          "Leave this blank for now",
                          "0 for 'cam0_image0.png', 1 for 'cam0image  0.bmp'",
                          "Number of cameras",
                          "0 for '0 deg', 1 for '90 deg', 2 for '180 deg', 3 for '270 deg'",
                          "The initial frame you want EventViewer to display",
                          "The trigger frame (defined by your DAQ)",
                          "Piezo variable name",
                          "Dytran variable name"]
    parameter_stringvars = []
    _max_width = max([len(p)+1 for p in parameter_labels])
    for n, parameter in enumerate(parameter_labels):
        p_frame = tk.Frame(master=canvas)
        p_frame.grid(row=len(path_labels)+n+3, column=1, sticky=tk.E)
        lab =tk.Label(master=p_frame, text=" "*(_max_width-len(parameter)-1)+parameter+":",
                      width= _max_width, font=DEFAULT_FONT)
        WidgetToolTip(lab, parameter_tooltips[n])
        lab.grid(row=0, column=0, sticky=tk.W)
        parameter_stringvars.append(tk.StringVar(master=p_frame, value=parameter_defaults[n]))
        tk.Entry(master=p_frame, width=50, textvariable=parameter_stringvars[n], font=DEFAULT_FONT)\
                .grid(row=0, column=1, sticky=tk.W)
    ## DISABLE RELATIVE IMAGE PATH (Note: .grid_slaves() returns in chronologically reversed grid order
    canvas.grid_slaves()[7].grid_slaves()[0].config(state=tk.DISABLED)
    canvas.grid_slaves()[7].grid_slaves()[1].config(state=tk.DISABLED)
    ## DISABLE RELATIVE IMAGE PATH (Note: .grid_slaves() returns in chronologically reversed grid order

    button_frame = tk.Frame(master=root)
    button_frame.grid(row=1, column=0, sticky=tk.E)
    cancel_button = tk.Button(master=button_frame, text="Cancel", font=DEFAULT_FONT, command=_on_cancel)
    cancel_button.grid(row=0, column=0, sticky=tk.W)
    ok_button = tk.Button(master=button_frame, text="OK", font=DEFAULT_FONT, command=_on_ok)
    ok_button.grid(row=0, column=1, sticky=tk.E)
    root.mainloop()

    file_paths = [path_stringvars[n].get() for n in range(len(path_labels))]
    return exit_code, file_paths


def run_intro():
    # Normally, I'd agree that global variables are bad and should be avoided,
    # however in the case of an exit_code I think the use is appropriate.
    global exit_code
    exit_code = 0
    def _on_ok():
        global exit_code
        exit_code = 1
        root.destroy()
        return

    def _on_cancel():
        root.destroy()
        return

    root = tk.Tk()
    root.eval("tk::PlaceWindow %s center" % root.winfo_pathname(root.winfo_id()))
    root.title(ROOT_TITLE)
    intro_h = 300; intro_w = 300

    canvas = tk.Canvas(master=root, height = intro_h, width = intro_w)
    root.resizable(False, False)


    _img = ImageTk.PhotoImage(Image.open("../Resources/placeholder.jpg"))
    canvas.create_image(0, 0, anchor=tk.NW, image=_img)
    canvas.pack()

    _text = "Welcome to EventViewer setup!\n"\
            "This tool will help you create\n"\
            "a config file for EventViewer.\n"\
            "     Press OK to continue"
    textx, texty = rel_to_abs(0.5, 0.5, intro_h, intro_w)
    canvas.create_text(textx, texty, text=_text, anchor=tk.CENTER, font=DEFAULT_FONT)

    button_frame = tk.Frame(master=canvas)
    button_frame.place(relx= 0.9, rely=0.9, anchor="se")

    cancel_button = tk.Button(master=button_frame,
                              text="Cancel", font=DEFAULT_FONT,
                              command=_on_cancel)
    cancel_button.grid(row=0, column=0, sticky=tk.E)
    ok_button = tk.Button(master=button_frame,
                          text="OK", font=DEFAULT_FONT,
                          command=_on_ok)
    ok_button.grid(row=0, column=1)

    root.mainloop()

    # It may be a little redundant to return exit_code since it is already in the global scope,
    # however the code is much clearer this way to anyone trying to debug it.
    return [exit_code, None]




def abs_to_rel(x, y, xmax, ymax):
    # Returns a scaled version of x and y between 0 and 1 where
    # xmax and ymax are the absolute sizes of the window.
    assert (x<=xmax and y<=ymax), \
        "x: {} and y: {} must be less than xmax: {} and ymax: {}".format(x, y, xmax, ymax)
    return [x/xmax, y/ymax]


def rel_to_abs(x, y, xmax, ymax):
    # Inverse of abs_to_rel.
    assert (0<=x<=1 and 0<=y<=1), \
        "x: {} and y: {} must each be between 0 and 1.".format(x, y)
    return [x*xmax, y*ymax]


if __name__ == "__main__":
    #exit_code, _ = run_setup()
    intro_test   = 0
    idk_test     = 1
    convert_test = 0

    if intro_test:
        run_intro()
    if idk_test:
        idk_yet()
    if convert_test:
        ask_and_convert_npy()


    #print("Setup completed with exit code: {}".format(exit_code))
