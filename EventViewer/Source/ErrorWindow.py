from tkinter.messagebox import showerror


class ErrorWindow:
    def __init__(self, master, message):
        showerror(master=master, title="Error!", message=message)
        return
