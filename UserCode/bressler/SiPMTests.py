
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
from scipy.fftpack import fft

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
        self.fig, self.ax= plt.subplots(nrows=2, ncols=1)
        self.runid = runid
        self.ev = ev
        self.areas=areas
        self.left = left
        self.right = right
        self.update_title()
        cid = plt.gcf().canvas.mpl_connect("key_press_event", self.key_press)

        return

    def plot_SiPM_trace(self):
        # plt.cla()
        self.ax[0].clear()
        self.ax[1].clear()
        xd = 1e9 * np.arange(pmt_data["traces"].shape[2]) * pmt_data["dt"][self.trigger.trig, 0]
        yd_0 = pmt_data["traces"][self.trigger.trig, 0, :] * pmt_data["v_scale"][self.trigger.trig, 0] + \
                  pmt_data["v_offset"][self.trigger.trig, 0]
        yd_1 = pmt_data["traces"][self.trigger.trig, 1, :] * pmt_data["v_scale"][self.trigger.trig, 1] + \
                  pmt_data["v_offset"][self.trigger.trig, 1]
        #fig,ax = plt.subplots(nrows=2, ncols=1)
        self.ax[0].plot(xd, yd_0)
        self.ax[0].plot(xd, yd_1)
        self.ax[0].set_xlabel("Time (ns)")
        self.ax[0].set_ylabel("Amp (a.u.)")

        yf0 = fft(yd_0)
        yf1 = fft(yd_1)
        xf = np.linspace(0.0, 1e9*len(xd)/(2*xd[-1]), len(xd)//2)
        self.ax[1].semilogy(xf, 2.0 / len(xd) * np.abs(yf0[0:len(xd) // 2]))
        self.ax[1].semilogy(xf, 2.0 / len(xd) * np.abs(yf1[0:len(xd) // 2]))
        self.ax[0].axvline(self.left, color="k", linewidth=2)
        self.ax[0].axvline(self.right, color="k", linewidth=2)
        dt = pmt_data["dt"][self.trigger.trig, 0]
        #plt.specgram(yd_1, NFFT=256, Fs=1e-9 /dt, )
        self.ax[1].set_xlabel("Frequency")
        self.ax[1].set_ylabel("Amplitude")
        self.update_title()
        plt.tight_layout()
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
    #raw_directory = r"C:\Users\John\Documents\SBC-18-data"
    raw_directory = "/bluearc/storage/SBC-18-data/"
    bias58vledON   = os.path.join(raw_directory, "20180906_5")
    bias586vledON  = os.path.join(raw_directory, "20180906_6")
    bias58vledOFF  = os.path.join(raw_directory, "20180906_7")
    bias586vledOFF = os.path.join(raw_directory, "20180906_8")
    bias59vledON   = os.path.join(raw_directory, "20180910_0")
    bias596vledON  = os.path.join(raw_directory, "20180910_1")
    bias587vledON  = os.path.join(raw_directory, "20180912_0")
    shortercables  = os.path.join(raw_directory, "20180914_0")
    nofan          = os.path.join(raw_directory, "20180918_0")
    nofan2         = os.path.join(raw_directory, "20180918_1")
    freddy         = os.path.join(raw_directory, "20180920_0")
    singlephotonpls= os.path.join(raw_directory, "20180919_0")
    v135           = os.path.join(raw_directory, "20180921_0")
    v136           = os.path.join(raw_directory, "20180921_1")
    v137           = os.path.join(raw_directory, "20180921_2")
    v138           = os.path.join(raw_directory, "20180921_3")
    v139           = os.path.join(raw_directory, "20180921_4")
    v140           = os.path.join(raw_directory, "20180921_5")
    v141           = os.path.join(raw_directory, "20180921_6")
    v142           = os.path.join(raw_directory, "20180921_7")
    v143           = os.path.join(raw_directory, "20180921_8")
    v144           = os.path.join(raw_directory, "20180921_9")
    new            = os.path.join(raw_directory, "20181205_3")
    labels = ["58V - LED ON",
              "58V - LED OFF",
              "58.6V - LED ON",
              "58.6V - LED OFF",
              "59V - LED ON",
              "59.6V - LED ON",
              "58.7V - LED ON (new)",
              "Shorter cables test",
              "no fan",
              "no fan 2",
              "freddy",
              "single photon maybe pls",
              "vhigh=1.35v",
              "vhigh=1.36v",
              "vhigh=1.37v",
              "vhigh=1.38v",
              "vhigh=1.39v",
              "vhigh=1.40v",
              "vhigh=1.41v",
              "vhigh=1.42v",
              "vhigh=1.43v",
              "vhigh=1.44v",
              "new"]

    var_array = [bias58vledON,
                 bias58vledOFF,
                 bias586vledON,
                 bias586vledOFF,
                 bias59vledON,
                 bias596vledON,
                 bias587vledON,
                 shortercables,
                 nofan,
                 nofan2,
                 freddy,
                 singlephotonpls,
                 v135,
                 v136,
                 v137,
                 v138,
                 v139,
                 v140,
                 v141,
                 v142,
                 v143,
                 v144,
                 new
                 ]

    colors = ["darkorange",
              "yellow",
              "cyan",
              "magenta",
              "green",
              "black",
              "blue",
              "red",
              "green",
              "darkorange",
              "magenta",
              "cyan",
              "red",
              "blue",
              "magenta",
              "green",
              "black",
              "yellow",
              "orange",
              "lime",
              "gray",
              "chocolate",
              "red"]

    active = np.array([0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0,
                       1], dtype=bool)

    labels = np.array(labels)[active]
    var_array = np.array(var_array)[active]
    colors = np.array(colors)[active]
    nbins = 2500
    #fig, ax = plt.subplots(1, 1)
    #plot_SiPM_trace(ax, SBCcode.get_event(active_event, 0, "PMTtraces")["PMTtraces"])

    # plt.ioff()
    areas = []
    max_times = []
    n_runs =  len(var_array)
    plt.ioff()
    left_lim = 800
    right_lim = 1400
    for run_ix in range(n_runs):
        sub_areas = []
        sub_max_times = []
        active_event = var_array[run_ix]
        events = SBCtools.BuildEventList(active_event)
        for ev in events:
            pmt_data = SBCcode.get_event(active_event, ev, "PMTtraces", max_file_size=1300)["PMTtraces"]
            n_triggers = pmt_data["traces"].shape[0]
            for trig in range(n_triggers):
                xd = 1e9 * np.arange(pmt_data["traces"].shape[2]) * pmt_data["dt"][trig, 0]
                yd_0 = pmt_data["traces"][trig, 0, :] * pmt_data["v_scale"][trig, 0] + \
                       pmt_data["v_offset"][trig, 0]
                yd_1 = pmt_data["traces"][trig, 1, :] * pmt_data["v_scale"][trig, 1] + \
                       pmt_data["v_offset"][trig, 1]
                good_indices = (xd < right_lim) & (xd > left_lim)
                avg = np.average(yd_0[:left_lim])
                #pmt_area = scipy.integrate.trapz(yd_1[good_indices]-avg,
                #                                 dx=1e9*pmt_data["dt"][trig, 1])
                # if np.any(yd_0 > 0.15):
                #     pmt_area = -1.
                # else:
                pmt_area = np.sum(yd_0[left_lim:right_lim]-avg)
                print(np.sum(yd_0))
                sub_areas.append(pmt_area)
                sub_max_times.append(1e9 * pmt_data["dt"][trig, 0] * np.argmax(pmt_data["traces"][trig, 0, :]))

            plotter = SiPMPlotter(pmt_data, runid=os.path.basename(active_event), ev=ev,
                                  areas=sub_areas,
                                  left=left_lim, right=right_lim)

            plotter.plot_SiPM_trace()
            plt.show()
        areas.append(sub_areas)
        max_times.append(sub_max_times)

    plt.hist(areas, bins=nbins, fill=False, color=colors, histtype="step", stacked=False, label=labels)
    plt.legend(loc="upper right")
    plt.suptitle("Areas")
    plt.yscale("log")

    # plt.figure()
    # plt.hist(max_times, bins=nbins, fill=False, color=colors, histtype="step", stacked=False)
    # plt.legend(labels)
    # plt.suptitle("Time of max")
    #plt.ion()
    plt.grid()
    plt.show()
    pass
