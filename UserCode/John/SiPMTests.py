
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate

import SBCcode
from SBCcode.Tools import SBCtools

class SiPMTrigger(object):
    def __init__(self, trig=0):
        self.trig=trig

class SiPMPlotter(object):
    def __init__(self, pmt_data, left=None, right=None,
                 default_trigger=0, runid=None, ev=None, areas=None):
        self.pmt_data = pmt_data
        self.default_trigger = default_trigger
        self.trigger = SiPMTrigger(trig=self.default_trigger)
        self.fig = plt.figure()
        self.runid = runid
        self.ev = ev
        self.areas=areas
        self.left = left
        self.right = right
        self.update_title()
        cid = plt.gcf().canvas.mpl_connect("key_press_event", self.key_press)

        return

    def plot_SiPM_trace(self):
        plt.cla()
        xd = 1e9 * np.arange(pmt_data["traces"].shape[2]) * pmt_data["dt"][self.trigger.trig, 0]
        yd_0 = pmt_data["traces"][self.trigger.trig, 0, :] * pmt_data["v_scale"][self.trigger.trig, 0] + \
                  pmt_data["v_offset"][self.trigger.trig, 0]
        yd_1 = pmt_data["traces"][self.trigger.trig, 1, :] * pmt_data["v_scale"][self.trigger.trig, 1] + \
                  pmt_data["v_offset"][self.trigger.trig, 1]
        plt.plot(xd, yd_0)
        plt.plot(xd, yd_1)
        plt.axvline(self.left, color="k", linewidth=2)
        plt.axvline(self.right, color="k", linewidth=2)
        dt = pmt_data["dt"][self.trigger.trig, 0]
        #plt.specgram(yd_1, NFFT=256, Fs=1e-9 /dt, )
        self.update_title()
        plt.draw()
        return

    def key_press(self, mpl_event):
        if mpl_event.key == "left":
            self.trigger.trig -= 1 if self.trigger.trig else -(self.pmt_data["traces"].shape[0]-1)
        elif mpl_event.key == "right":
            self.trigger.trig += 1 if self.trigger.trig < self.pmt_data["traces"].shape[0] else -(self.pmt_data["traces"].shape[0]-1)
        else:
            return
        self.plot_SiPM_trace()
        return

    def update_title(self):
        plt.suptitle("Runid: {} || Event: {} || Trigger {}\nUse Left and Right arrows to navigate.\nArea={}"\
                     .format(self.runid, self.ev, self.trigger.trig, self.areas[self.trigger.trig]))

if __name__ == "__main__":
    raw_directory = r"C:\Users\John\Documents\SBC-18-data"
    bias58vledON   = os.path.join(raw_directory, "20180906_5")
    bias586vledON  = os.path.join(raw_directory, "20180906_6")
    bias58vledOFF  = os.path.join(raw_directory, "20180906_7")
    bias586vledOFF = os.path.join(raw_directory, "20180906_8")
    labels = ["58V - LED ON",
              "58V - LED OFF",
              "58.6V - LED ON",
              "58.6V - LED OFF"]
    var_array = [bias58vledON,
                 bias58vledOFF,
                 bias586vledON,
                 bias586vledOFF]
    colors = ["darkorange",
              "red",
              "cyan",
              "magenta"]
    active = np.array([1,0,1,0], dtype=bool)
    labels = np.array(labels)[active]
    var_array = np.array(var_array)[active]
    colors = np.array(colors)[active]

    #fig, ax = plt.subplots(1, 1)
    #plot_SiPM_trace(ax, SBCcode.get_event(active_event, 0, "PMTtraces")["PMTtraces"])

    # plt.ioff()
    areas = []
    max_times = []
    n_runs =  len(var_array)
    plt.ioff()
    left_lim = 700
    right_lim = 1300
    for run_ix in range(n_runs):
        sub_areas = []
        sub_max_times = []
        active_event = var_array[run_ix]
        events = SBCtools.BuildEventList(active_event)
        for ev in events:
            pmt_data = SBCcode.get_event(active_event, ev, "PMTtraces")["PMTtraces"]
            n_triggers = pmt_data["traces"].shape[0]
            for trig in range(n_triggers):
                xd = 1e9 * np.arange(pmt_data["traces"].shape[2]) * pmt_data["dt"][trig, 0]
                yd_0 = pmt_data["traces"][trig, 0, :] * pmt_data["v_scale"][trig, 0] + \
                       pmt_data["v_offset"][trig, 0]
                yd_1 = pmt_data["traces"][trig, 1, :] * pmt_data["v_scale"][trig, 1] + \
                       pmt_data["v_offset"][trig, 1]
                good_indices = (xd < right_lim) & (xd > left_lim)

                pmt_area = scipy.integrate.trapz(yd_0[good_indices],
                                                 dx=1e9*pmt_data["dt"][trig, 1])

                sub_areas.append(pmt_area)
                sub_max_times.append(1e9 * pmt_data["dt"][trig, 0] * np.argmax(pmt_data["traces"][trig, 0, :]))
                pass
            plotter = SiPMPlotter(pmt_data, runid=os.path.basename(active_event), ev=ev,
                                  areas=sub_areas,
                                  left=left_lim, right=right_lim)
            plotter.plot_SiPM_trace()

            plt.show()
        areas.append(sub_areas)
        max_times.append(sub_max_times)

    plt.hist(areas, bins=3000, fill=False, color=colors, histtype="step", stacked=False)
    plt.legend(labels)
    plt.suptitle("Areas")
    plt.figure()
    plt.hist(max_times, bins=3000, fill=False, color=colors, histtype="step", stacked=False)
    plt.legend(labels)
    plt.suptitle("Time of max")
    plt.show()
    pass
