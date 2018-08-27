## John Gresl



import numpy as np
import scipy
import re
import sys
import time
import scipy.signal
import fcntl
import scipy.optimize
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
try:
    from PyQt4 import QtCore
    PyQt4import = True
except:
    PyQt4import = False
import os
from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE



class mplActions(object):
    def __init__(self, timebase, piezo, event_id, event, out_file_obj, drawing_ax):
        self.timebase = timebase
        self.piezo = piezo
        self.out_file_obj = out_file_obj
        self.event_id = event_id
        self.event = event
        self.drawing_ax = drawing_ax
        self.template_str = "{}, {}, {}\n"
        self.active_t0 = None
        self.full_stop = False
        return

    def on_click(self, mpl_event):
        if not mpl_event.inaxes == self.drawing_ax:
            return
        if not mpl_event.dblclick:
            return
        cur_xlim = self.drawing_ax.get_xlim()
        cur_ylim = self.drawing_ax.get_ylim()
        self.drawing_ax.cla()
        self.drawing_ax.plot(self.timebase, self.piezo)
        self.drawing_ax.axvline(x=mpl_event.xdata, color="k", linewidth=2)
        self.drawing_ax.set_xlim(cur_xlim)
        self.drawing_ax.set_ylim(cur_ylim)
        plt.draw()
        self.active_t0 = mpl_event.xdata
        return

    def confirm(self, click_ev):
        if self.active_t0 is None:
            print("Click the plot to select a t0 first!")
            return
        else:
            print("t0 = {}".format(self.active_t0))
            self.out_file_obj.write(self.template_str.format(self.event_id,
                                                             self.event,
                                                             self.active_t0))
        plt.close()
        return

    def skip(self, click_ev):
        print("Skipping. Not writing anything")
        plt.close()
        return

    def stop(self, click_ev):
        print("Closing! Bye bye.")
        self.full_stop = True
        plt.close()
        return



def run(data_directory, output_file, data_set = "SBC-17", start=0, n_events=5):

    ## Want to go through events and show the raw piezo trace, then click to select a t0.
    ## Display the t0, ask for confirmation, and then save to an output file.
    full_stop = False

    with open(output_file, "w") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB) # <-- Locks the file until f.__exit__()
        except BlockingIOError:
            print("Critical: File is currently locked by another user/process. Aborting.")
            return
        except:
            print("Unsupported on non-Linux/Mac OS's. Aborting.")
            return
        evs_processed = 0
        all_ids = [i for i in os.listdir(data_directory) if
                   os.path.isdir(os.path.join(data_directory, i)) and i != "scratch"]
        n_skip = 0
        _counter = 0
        for ev_id in all_ids:
            if data_set == "SBC-17":
                f_template = re.compile("^[0-9]{8}_[0-9]+")  ## Until the year 10,000 this should work. . .
                f_match = re.match(f_template, ev_id)
                if f_match:
                    d = int(ev_id.split("_")[0])
                    if d < 20170630:
                        continue  # Skips all data that happened before June 30 2017.
                all_evs = [i for i in os.listdir(os.path.join(data_directory, ev_id)) if
                           os.path.isdir(os.path.join(data_directory, ev_id, i))]
                all_evs = sorted(all_evs, key=int)
                for ev in all_evs:
                    if start > _counter:
                        _counter += 1
                        continue
                    if full_stop:
                        # f.write(_last_out_str)
                        print("Processed {} events ({}-{})".format(evs_processed - 1, start, start + evs_processed - 1))
                        return
                    if evs_processed >= n_events:
                        # f.write(_last_out_str)
                        print("Processed {} events ({}-{})".format(evs_processed, start, start + evs_processed))
                        return


                    else:
                        ## We grab the event and extract the Piezo traces
                        print(os.path.join(RAW_DIRECTORY, ev_id, ev), end=" >>>>> ")
                        ev_data = GE(os.path.join(RAW_DIRECTORY, ev_id), ev, "fastDAQ")
                        try:
                            piezo2 = ev_data["fastDAQ"]["Piezo2"]
                        except Exception as e:
                            print("********************Exception! Fix this.********************")
                            continue

                        ## We grab the timebase from the event data and find the mean time interval
                        timebase = ev_data["fastDAQ"]["time"]
                        xextent = (min(timebase), max(timebase))
                        dt = np.mean(np.diff(timebase))
                        plt.plot(timebase, piezo2)
                        drawing_ax = plt.gca()
                        plt.ioff()
                        actions = mplActions(timebase, piezo2, ev_id, ev, f, drawing_ax)
                        cid = plt.gcf().canvas.mpl_connect('button_press_event', actions.on_click)

                        g_fit = plt.axes([0.543, 0.05, 0.357, 0.15])
                        g_but = Button(g_fit, "Confirm", hovercolor="green")
                        g_but.label.set_fontsize(45)
                        g_but.on_clicked(actions.confirm)

                        b_fit = plt.axes([0.125, 0.05, 0.359, 0.15])
                        b_but = Button(b_fit, "Skip", hovercolor="red")
                        b_but.label.set_fontsize(45)
                        b_but.on_clicked(actions.skip)

                        close_fit = plt.axes([0.91, 0.05, .08, 0.15])
                        close_button = Button(close_fit, "Stop\nall", hovercolor="red")
                        close_button.label.set_fontsize(14)
                        close_button.on_clicked(actions.stop)

                        plt.show()
                        full_stop = actions.full_stop
                        evs_processed += 1
    return




if __name__ == "__main__":
    RAW_DIRECTORY = "/bluearc/storage/SBC-17-data/"
    start = 0
    n_events = 3
    output_file = "{}_{}-{}.txt".format("johns_t0_handscan",start, start+n_events)
    run(RAW_DIRECTORY, output_file, start=start, n_events=n_events)

    pass