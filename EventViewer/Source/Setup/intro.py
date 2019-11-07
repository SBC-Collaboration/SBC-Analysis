import tkinter as tk
from EventViewer.Source.GUIDefaults import *

class EVIntro(object):
    def __init__(self):
        self.ok_clicked = False # After destruction of tk objects, the outside API should check this parameter
                                # to see if the user pressed ok or not and continue accordingly.
        self.root = tk.Tk()
        self.root.eval("tk::PlaceWindow %s center" % \
                       self.root.winfo_pathname(self.root.winfo_id()))
        self.root.title("EventViewer Setup")
        self.canvas = tk.Canvas(master=self.root,
                                background=BACKGROUND)
        self.canvas.pack()
        #self.root.resizable(False, False)
        text = "\n\nWelcome to EventViewer setup!\n" \
               "This tool will help you create\n" \
               "a config file for EventViewer.\n" \
               "     Press OK to continue.     "
        text_tk = tk.Text(master=self.canvas, font=(FONT_NAME, FONT_SIZE),
                          background=BACKGROUND, height=6, width=30)
        text_tk.tag_configure("center", justify="center")
        text_tk.insert("1.0", text)
        text_tk.tag_add("center", "1.0", "end")
        text_tk.pack()


        #self.root.maxsize(320, 300)
        #self.root.minsize(320, 300)
        self.ok_button = tk.Button(master=self.canvas, text="Ok",
                                   font=(FONT_NAME, FONT_SIZE), background=BACKGROUND,
                                   command=self.ok)
        self.ok_button.pack()
        self.cancel_button = tk.Button(master=self.canvas, text="Cancel",
                                       font=(FONT_NAME, FONT_SIZE), background=BACKGROUND,
                                       command=self.cancel)
        self.cancel_button.pack()
        self.root.resizable(False, False)

        self.root.mainloop()
        return

    def cancel(self):
        self.root.destroy()
        return

    def ok(self):
        self.root.destroy()
        self.ok_clicked = True
        return


def main():
    test = EVIntro()
    if not test.ok_clicked:
        return




if __name__ == "__main__":
    main()

