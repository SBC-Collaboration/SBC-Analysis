from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk


class CustomToolbar(NavigationToolbar2Tk):
    # Creates a custom toolbar that only uses the buttons we want. Ideally,
    # this removes the forward, backward, and configure subplots button.
    toolitems = [t for t in NavigationToolbar2Tk.toolitems] #3if \
        #t[0] in ("Home", "Pan", "Zoom", "Save")]
