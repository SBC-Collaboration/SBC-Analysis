#! /coupp/app/home/coupp/anaconda3/bin/python3
## John Gresl

##



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


matplotlib.use("Qt4Agg")


def extend_window(w, r):
    # Extends a "window" by a ratio "r".
    mp = 0.5*(w[1]+w[0])  ## Midpoint
    new_len = (w[1]-w[0])*(1+r) ## Length of new window
    return [mp-new_len/2, mp+new_len/2]


class FitCheckButtonActions(object):
    def __init__(self, event_id, ev, tau, out_file_obj):
        self.event_id = event_id
        self.ev = ev
        self.tau = tau
        self.out_file_obj = out_file_obj
        self.template_str = "{}, {}, {}\n"
        self.full_stop = False
        return

    def good_fit(self, click_ev):
        print("Good fit!")
        self.out_file_obj.write(self.template_str.format(self.event_id,
                                                         self.ev,
                                                         self.tau))
        plt.close()
        return

    def bad_fit(self, click_ev):
        print("bad fit...")
        self.out_file_obj.write(self.template_str.format(self.event_id,
                                                         self.ev,
                                                         -1))
        plt.close()
        return

    def stop(self, click_ev):
        print("Closing! Bye bye.")
        self.full_stop = True
        plt.close()
        return


def calculate_sums(Pxx, good_indices, n):
    out = []
    for j in range(n):
        out.append(sum(Pxx[good_indices, j]))
    return out


def rescale_window(w1, w2):
    # Returns a scaled version of w2 so that the endpoints of w2 match those of w1
    y1, y2 = min(w1), max(w1)
    x1, x2 = min(w2), max(w2)
    a = (y1-y2)/(x1-x2)
    b = (x1*y2-x2*y1)/(x1-x2)
    return a*w2+b


def _bandpass(data, lower=None, upper=None):
    if lower is None and upper is None:
        return data
    if lower is None:
        return np.where([x <=  upper for x in data])
    if upper is None:
        return np.where([x >= lower for x in data])
    return np.where([lower <= x <= upper for x in data])


def _spectrum_sums(spectrum, good_indices, N):
    out = []
    for n in range(N):
        out.append(sum(spectrum[good_indices[0], n]))
    return out


def remove_ticks(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    return


def is_int(n):
    try:
        int(n)
        return True
    except:
        return False


def t0_finder(spectrum_sums, timebase):
    ## A first attempt might be to find the maximum in the spectrum and use that as the t0
    ind = np.argmax(spectrum_sums)
    return timebase[ind]


def exponential_fit(x, y, t0=None, lookahead=0.02):
    if t0 is None:
        t0 = t0_finder(y, x)
    left_ind = np.argmin(np.abs(x-t0))
    right_ind = np.argmin(np.abs(x-(t0+lookahead)))
    fit = scipy.optimize.curve_fit(lambda t, a, b: a*np.exp(-b*(t-t0)),
                                   x[left_ind:right_ind], y[left_ind:right_ind])
    xout = x[left_ind:right_ind]
    yout = fit[0][0]*np.exp(-fit[0][1]*(xout-t0))
    return (xout, yout, fit[0][0], fit[0][1])


def run(data_directory, output_file, data_set = "SBC-17", start=0, n_events=5):
    ## Inputs:
    ##      data_directory: Directory where the raw data is
    ##      output_file:    Location to store output binary file.
    ##      data_set:       String signifying which dataset we're using so we can filter out known 'bad-data'...
    ##      start:          Starts on event number 'start
    ##      n_events:       Number of events to process
    LOWER = 15000
    UPPER = 40000
    ## LOWER and UPPER are the frequency limits we place on the Spectogram
    #if os.path.isfile(output_file):
    #    already_scanned = RB(output_file)
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


        already_scanned = {"run_id": [],
                           "tau": []}
        evs_processed = 0
        all_ids = [i for i in os.listdir(data_directory) if os.path.isdir(os.path.join(data_directory, i)) and i != "scratch"]
        n_skip = 0
        _counter = 0
        for ev_id in all_ids:
            if data_set == "SBC-17":
                f_template = re.compile("^[0-9]{8}_[0-9]+") ## Until the year 10,000 this should work. . .
                f_match = re.match(f_template, ev_id)
                if f_match:
                    d = int(ev_id.split("_")[0])
                    if d < 20170630:
                        continue # Skips all data that happened before June 30 2017.
                all_evs = [i for i in os.listdir(os.path.join(data_directory, ev_id)) if os.path.isdir(os.path.join(data_directory, ev_id, i))]
                all_evs = sorted(all_evs, key=int)
                for ev in all_evs:
                    if start > _counter:
                        _counter += 1
                        continue
                    if full_stop:
                        #f.write(_last_out_str)
                        print("Processed {} events ({}-{})".format(evs_processed-1, start, start + evs_processed-1))
                        return
                    if  evs_processed >= n_events:
                        #f.write(_last_out_str)
                        print("Processed {} events ({}-{})".format(evs_processed, start, start + evs_processed))
                        return
                    if (ev_id, ev) in already_scanned["run_id"]:
                        n_skip += 1
                        print("Skipping", (ev_id, ev), ". Previously scanned.")

                    else:
                        ## We grab the event and extract the Piezo traces
                        print(os.path.join(RAW_DIRECTORY, ev_id, ev), end=" >>>>> ")
                        ev_data = GE(os.path.join(RAW_DIRECTORY, ev_id), ev, "fastDAQ")
                        try:
                            piezo1 = ev_data["fastDAQ"]["Piezo1"]
                            piezo2 = ev_data["fastDAQ"]["Piezo2"]
                        except Exception as e:
                            print("********************Exception! Fix this.********************")
                            continue

                        ## We grab the timebase from the event data and find the mean time interval
                        timebase = ev_data["fastDAQ"]["time"]
                        xextent = (min(timebase), max(timebase))
                        dt = np.mean(np.diff(timebase))

                        ## We start to create our plots and name them for future use
                        fig, axarr = plt.subplots(4,2, figsize=(15,9), sharex=True)
                        spectrum1_ax = axarr[0][0]
                        spectrum2_ax = axarr[0][1]
                        raw1_ax      = axarr[1][0]
                        raw2_ax      = axarr[1][1]
                        filtered1_ax = axarr[2][0]
                        filtered2_ax = axarr[2][1]
                        fit1_ax      = axarr[3][0]
                        fit2_ax      = axarr[3][1]
                        #conv1_ax     = axarr[4][0]
                        #conv2_ax     = axarr[4][1]

                        statistics1_axes = axarr[-1][0]
                        statistics2_axes = axarr[-1][1]


                        ## Create the specgrams and extract useful information such as the spectrum (Pxx), the frequency bins, the time
                        ## interval bins, and the matplotlib image object
                        spectrum1, freqs1, bins1, im1 = spectrum1_ax.specgram(piezo1, Fs=1./dt, NFFT=256, xextent=xextent, mode="magnitude")
                        spectrum2, freqs2, bins2, im2 = spectrum2_ax.specgram(piezo2, Fs=1./dt, NFFT=256, xextent=xextent, mode="magnitude")
                        spectrum1_ax.set_title("Piezo1 -- Disabled since we only use 2 for now")
                        spectrum2_ax.set_title("Piezo2")
                        spectrum1_ax.set_ylim(extend_window([LOWER, UPPER], 0.05))
                        spectrum2_ax.set_ylim(extend_window([LOWER, UPPER], 0.05))

                        ## Plot the raw signal and set some labels/other aesthetic choices
                        raw1_ax.plot(timebase, piezo1)
                        raw2_ax.plot(timebase, piezo2)

                        ## Now we start our signal filtering!
                        ## First step is to find the indices of all the frequencies we are interested in. These are
                        ## the frequencies between LOWER and UPPER defined above.
                        good_indices1 = _bandpass(freqs1, LOWER, UPPER)
                        good_indices2 = _bandpass(freqs2, LOWER, UPPER)

                        ## Now we go through and calculate the sums of the fft power spectrum for each time, only counting the
                        ## frequencies we are interested in
                        N1 = len(bins1)
                        N2 = len(bins2) # In theory N1=N2 because of the sample rate of the DAQ but rather safe than sorry.
                        spectrum1_sums = _spectrum_sums(spectrum1, good_indices1, N1)
                        spectrum2_sums = _spectrum_sums(spectrum2, good_indices2, N2)

                        ## Now we apply a median filter to the data to attempt to remove some noise.
                        spectrum1_sums = scipy.signal.medfilt(spectrum1_sums)
                        spectrum2_sums = scipy.signal.medfilt(spectrum2_sums)

                        ## Rescale the time window so we can plot it appropriately.
                        rescaled_time1 = rescale_window(xextent, bins1)
                        rescaled_time2 = rescale_window(xextent, bins2)

                        ## Plot the spectrum sums.
                        filtered1_ax.plot(rescaled_time1, spectrum1_sums)
                        filtered2_ax.plot(rescaled_time2, spectrum2_sums)

                        ## Extract t0 from previous analysis
                        reco_data = RB(os.path.join(RECO_DIRECTORY, "{}/AcousticAnalysis_{}.bin".format(ev_id, ev_id)))
                        old_t0 = sum(reco_data["bubble_t0"][int(ev)])/2

                        ## Attempt to independently determine t0 from the filtered spectrum sum
                        new_t0_1 = t0_finder(spectrum1_sums, rescaled_time1)
                        new_t0_2 = t0_finder(spectrum2_sums, rescaled_time2)

                        ## Exponential fit
                        # scipy.optimize.curve_fit(lambda t,a,b: a*numpy.exp(b*t),  x,  y)
                        fitx1, fity1, a1, b1 = exponential_fit(rescaled_time1, spectrum1_sums, new_t0_1)
                        fitx2, fity2, a2, b2 = exponential_fit(rescaled_time2, spectrum2_sums, new_t0_2)

                        filtered1_ax.axvspan(fitx1[0], fitx1[-1], facecolor="g", alpha=0.5)
                        filtered2_ax.axvspan(fitx2[0], fitx2[-1], facecolor="g", alpha=0.5)
                        filtered1_ax.plot(fitx1, fity1, "r-", linewidth=3)
                        filtered2_ax.plot(fitx2, fity2, "r-", linewidth=3)

                        fit1_ax.plot(fitx1, fity1,"r-", linewidth=3)
                        fit2_ax.plot(fitx2, fity2, "r-", linewidth=3)
                        fit1_ax.plot(rescaled_time1, spectrum1_sums)
                        fit2_ax.plot(rescaled_time2, spectrum2_sums)
                        fit1_ax.set_xlim([fitx1[0], fitx1[-1]])
                        fit2_ax.set_xlim([fitx2[0], fitx2[-1]])


                        ## Plot previous_t0 on all of the axes
                        for ax in axarr.flatten()[:-2]:
                            #for ax in inner:
                            #ax.set_xlim(extend_window([old_t0*(1+0.03), old_t0*(1-0.03)], 0.1))
                            ax.set_xlim(xextent)
                            ax.axvline(x=old_t0, color="k", linewidth=1)


                        ## Show some statistics in the final plot
                        statistics1_axes.text(0.5, 0.5,
                                              "Old t0: {}\nNew t0: {}\n% Diff: {:.5f}%\nTau: {}".\
                                              format(old_t0,
                                                     new_t0_1,
                                                     100*abs((old_t0 - new_t0_1) / old_t0),
                                                     1. / b1),
                                              fontsize=12, transform=statistics1_axes.transAxes,
                                              horizontalalignment="center", verticalalignment="center")
                        remove_ticks(spectrum2_ax)
                        remove_ticks(raw1_ax)
                        remove_ticks(raw2_ax)
                        remove_ticks(filtered1_ax)
                        remove_ticks(filtered2_ax)
                        remove_ticks(fit1_ax)
                        remove_ticks(fit2_ax)
                        #remove_ticks(conv1_ax)
                        #remove_ticks(conv2_ax)
                        remove_ticks(statistics1_axes)
                        remove_ticks(statistics2_axes)

                        plt.subplots_adjust(bottom=0.2, hspace=0, wspace=0)


                        statistics2_axes.text(0.5, 0.5,
                                              "Old t0: {}\nNew t0: {}\n% Diff: {:.5f}%\nTau: {}". \
                                              format(old_t0,
                                                     new_t0_2,
                                                     100 * abs((old_t0 - new_t0_2) / old_t0),
                                                     1. / b2),
                                              fontsize=12, transform=statistics2_axes.transAxes,
                                              horizontalalignment="center", verticalalignment="center")
                        actions = FitCheckButtonActions(event_id=ev_id, ev=ev, tau=1./b2, out_file_obj=f)
                        g_fit = plt.axes([0.543, 0.05, 0.357, 0.15])
                        g_but = Button(g_fit, "Piezo2\nGood Fit", hovercolor="green")
                        g_but.label.set_fontsize(45)
                        g_but.on_clicked(actions.good_fit)

                        b_fit = plt.axes([0.125,0.05,0.359,0.15])
                        b_but = Button(b_fit, "Piezo2\nBad Fit", hovercolor="red")
                        b_but.label.set_fontsize(45)
                        b_but.on_clicked(actions.bad_fit)

                        close_fit = plt.axes([0.91, 0.05, .08, 0.15])
                        close_button = Button(close_fit, "Stop\nall", hovercolor="red")
                        close_button.label.set_fontsize(14)
                        close_button.on_clicked(actions.stop)

                        ## Title the figure
                        fig.suptitle("Run: {}  ---  Event: {}".format(ev_id, ev))
                        win = plt.gcf().canvas.manager.window
                        if PyQt4import:
                            win.setWindowFlags(win.windowFlags() | QtCore.Qt.CustomizeWindowHint)
                            win.setWindowFlags(win.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)
                        plt.show()
                        full_stop = actions.full_stop
                        evs_processed += 1
        return # <-- This return shouldn't be reached since our exit condition is above




if __name__ == "__main__":

    RAW_DIRECTORY = "/bluearc/storage/SBC-17-data/"
    RECO_DIRECTORY = "/pnfs/coupp/persistent/grid_output/SBC-17/output/"
    OUTPUT_FILE = "SAMPLE.TXT"
    #run(RAW_DIRECTORY, OUTPUT_FILE, start=1520, n_events=80)
    run(RAW_DIRECTORY, OUTPUT_FILE, start=0, n_events=50)
## Jump forward about 10 ms, then look backward, if you see a cluster (another PMT pulse within 1 ms), discount that
## pmt pulse