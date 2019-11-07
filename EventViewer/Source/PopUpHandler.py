

import logging
from tkinter import messagebox


class PopUpHandler(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
        self.setFormatter(logging.Formatter(fmt='%(message)s'))
    def emit(self, message):
        messagebox.showerror('error', self.format(message))
