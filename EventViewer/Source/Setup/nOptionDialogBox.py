import tkinter as tk
import math

from EventViewer.Source.GUIDefaults import *

class nOptionDialogBox(object):
    def __init__(self, prompt=None, options=None, commands=None, master=None, ncols=4):

        self.prompt = "Please select an option." if prompt is None else prompt
        self.options = ["Yes", "No"] if options is None else options

        if master is None:
            self.root = tk.Tk()
            self.root.resizable(False, False)
        else:
            self.root = tk.Toplevel(master=master)
        self.commands = [root.destroy, root.destroy] if commands is None else commands
        if len(self.options) != len(self.commands):
            raise ValueError("There must be one command for each option and vice versa."+
                             "Got n_options={} and n_commands={}".format(len(self.options), len(self.commands)))
        self.selected_ix = None
        n = len(self.options)
        self.ncols = ncols if n >= ncols else n
        self.prompt_label = tk.Label(master=self.root, text=self.prompt, font=HEADER_FONT)
        self.prompt_label.grid(row=0, column=0, sticky=tk.NSEW, pady = 25)

        cur_col = 0
        cur_row = 0
        self.option_frame = tk.Frame(master=self.root)
        self.option_frame.grid(row=1, column=0, sticky=tk.EW)
        for opt in self.options:
            button = tk.Button(master=self.option_frame, text=opt, font=DEFAULT_FONT,
                               command=lambda col=cur_col, row=cur_row: self.on_select(col, row))
            button.grid(row=cur_row, column=cur_col, sticky=tk.EW, padx=10, pady=10)
            cur_col += 1
            if cur_col >= ncols:
                cur_row += 1
                cur_col = 0
        self.root.grab_set()
        return


    def on_select(self, col, row):
        self.selected_ix = self.ncols*row+col
        self.commands[self.selected_ix](self.selected_ix)
        self.root.grab_release()
        self.root.destroy()
        return


if __name__ == "__main__":

    def test(obj):
        print(obj)
        return


    def click_me(obj):
        options = list("ABC")
        funcs = []
        test = nOptionDialogBox(prompt="Pick your favorite letter", options=options, master=obj, commands=funcs)
        return

    root = tk.Tk()
    b = tk.Button(master=root, text="Click me", command=lambda obj=root: click_me(obj))
    b.pack()
    root.mainloop()






