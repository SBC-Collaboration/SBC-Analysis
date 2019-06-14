import os
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, freqz

#For Gaussian Fitting
from scipy.optimize import curve_fit
from numpy import exp, linspace, random

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

#Gaussian Function
def gaussian(x, amp, cen, wid):
    return amp * exp(-(x-cen)**2 / (wid))

#Sum Gaussian Function
def sum_gaussian(x, amp1, cen1, wid1, amp2, cen2, wid2):
    return (amp1 * exp(-(x-cen1)**2 / (wid1)))+ (amp2 * exp(-(x-cen2)**2 / (wid2)));

#Butterworth filter
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


if __name__ == "__main__":
    #raw_directory = r"C:\Users\John\Documents\SBC-18-data"
    # raw_directory = "/bluearc/storage/SBC-18-data/"
    # bias58vledON   = os.path.join(raw_directory, "20180906_5")
    # bias586vledON  = os.path.join(raw_directory, "20180906_6")
    # bias58vledOFF  = os.path.join(raw_directory, "20180906_7")
    # bias586vledOFF = os.path.join(raw_directory, "20180906_8")
    # bias59vledON   = os.path.join(raw_directory, "20180910_0")
    # bias596vledON  = os.path.join(raw_directory, "20180910_1")
    # bias587vledON  = os.path.join(raw_directory, "20180912_0")
    # shortercables  = os.path.join(raw_directory, "20180914_0")
    # nofan          = os.path.join(raw_directory, "20180918_0")
    # nofan2         = os.path.join(raw_directory, "20180918_1")
    # freddy         = os.path.join(raw_directory, "20180920_0")
    # singlephotonpls= os.path.join(raw_directory, "20180919_0")
    # v135           = os.path.join(raw_directory, "20180921_0")
    # v136           = os.path.join(raw_directory, "20180921_1")
    # v137           = os.path.join(raw_directory, "20180921_2")
    # v138           = os.path.join(raw_directory, "20180921_3")
    # v139           = os.path.join(raw_directory, "20180921_4")
    # v140           = os.path.join(raw_directory, "20180921_5")
    # v141           = os.path.join(raw_directory, "20180921_6")
    # v142           = os.path.join(raw_directory, "20180921_7")
    # v143           = os.path.join(raw_directory, "20180921_8")
    # v144           = os.path.join(raw_directory, "20180921_9")
    # new            = os.path.join(raw_directory, "20181205_3")
    # test1          = "20190125_1"
    # test2          = "20190222_0"

    raw_directory = "Data/"

    #Original runs
    zero = os.path.join(raw_directory,"20190516/20190516_0")
    one = os.path.join(raw_directory,"20190516/20190516_1")
    two = os.path.join(raw_directory,"20190516/20190516_2")
    three = os.path.join(raw_directory,"20190516/20190516_3")
    four = os.path.join(raw_directory,"20190516/20190516_4")
    five = os.path.join(raw_directory,"20190516/20190516_5")
    six = os.path.join(raw_directory,"20190516/20190516_6")
    seven = os.path.join(raw_directory,"20190516/20190516_7")
    eight = os.path.join(raw_directory,"20190516/20190516_8")
    nine = os.path.join(raw_directory,"20190516/20190516_9")
    ten = os.path.join(raw_directory,"20190516/20190516_10")
    eleven = os.path.join(raw_directory,"20190516/20190516_11")
    twelve = os.path.join(raw_directory,"20190516/20190516_12")
    thirteen = os.path.join(raw_directory,"20190516/20190516_13")
    fourteen = os.path.join(raw_directory,"20190516/20190516_14")
    fifteen = os.path.join(raw_directory,"20190516/20190516_15")
    sixteen = os.path.join(raw_directory,"20190516/20190516_16")
    seventeen = os.path.join(raw_directory,"20190516/20190516_17")
    eighteen = os.path.join(raw_directory,"20190516/20190516_18")
    nineteen = os.path.join(raw_directory,"20190516/20190516_19")
    twenty = os.path.join(raw_directory,"20190516/20190516_20")
    twentyone = os.path.join(raw_directory,"20190516/20190516_21")
    twentytwo = os.path.join(raw_directory,"20190516/20190516_22")
    twentythree = os.path.join(raw_directory,"20190516/20190516_23")
    twentyfour = os.path.join(raw_directory,"20190516/20190516_24")
    twentyfive = os.path.join(raw_directory,"20190516/20190516_25")
    twentysix = os.path.join(raw_directory,"20190516/20190516_26")
    twentyseven = os.path.join(raw_directory,"20190516/20190516_27")
    twentyeight = os.path.join(raw_directory,"20190516/20190516_28")

    #runs with reduced noise hardware
    clean17 = os.path.join(raw_directory, "20190523/20190523_17")
    clean18 = os.path.join(raw_directory, "20190523/20190523_18")
    clean19 = os.path.join(raw_directory, "20190523/20190523_19")

    #-166 temp runs
    Low_15 = os.path.join(raw_directory, "20190528/20190528_15")
    Low_16 = os.path.join(raw_directory, "20190528/20190528_16")
    Low_17 = os.path.join(raw_directory, "20190528/20190528_17")
    Low_18 = os.path.join(raw_directory, "20190528/20190528_18")
    Low_19 = os.path.join(raw_directory, "20190528/20190528_19")
    Low_20 = os.path.join(raw_directory, "20190528/20190528_20")
    Low_21 = os.path.join(raw_directory, "20190528/20190528_21")
    Low_22 = os.path.join(raw_directory, "20190528/20190528_22")
    Low_23 = os.path.join(raw_directory, "20190528/20190528_23")
    Low_24 = os.path.join(raw_directory, "20190528/20190528_24")
    Low_25 = os.path.join(raw_directory, "20190528/20190528_25")
    Low_26 = os.path.join(raw_directory, "20190528/20190528_26")
    Low_27 = os.path.join(raw_directory, "20190528/20190528_27")
    Low_28 = os.path.join(raw_directory, "20190528/20190528_28")

    #Dark Count Runs
    Dark_0 = os.path.join(raw_directory, "20190530/20190530_0")
    Dark_3 = os.path.join(raw_directory, "20190530/20190530_3")
    Dark_5 = os.path.join(raw_directory, "20190530/20190530_5")
    Dark_6 = os.path.join(raw_directory, "20190530/20190530_6")


    #
    # labels = ["58V - LED ON",
    #           "58V - LED OFF",
    #           "58.6V - LED ON",
    #           "58.6V - LED OFF",
    #           "59V - LED ON",
    #           "59.6V - LED ON",
    #           "58.7V - LED ON (new)",
    #           "Shorter cables test",
    #           "no fan",
    #           "no fan 2",
    #           "freddy",
    #           "single photon maybe pls",
    #           "vhigh=1.35v",
    #           "vhigh=1.36v",
    #           "vhigh=1.37v",
    #           "vhigh=1.38v",
    #           "vhigh=1.39v",
    #           "vhigh=1.40v",
    #           "vhigh=1.41v",
    #           "vhigh=1.42v",
    #           "vhigh=1.43v",
    #           "vhigh=1.44v",
    #           "new",
    #           "test1",
    #           "test2"]

    labels = ["zero",
            "one",
            "two",
            "three",
            "four",
            "five",
            "six",
            "seven",
            "eight",
            "nine",
            "ten",
            "eleven",
            "twelve",
            "thirteen",
            "fourteen",
            "fifteen",
            "sixteen",
            "seventeen",
            "eighteen",
            "nineteen",
            "twenty",
            "twentyone",
            "twentytwo",
            "twentythree",
            "twentyfour",
            "twentyfive",
            "twentysix",
            "twentyseven",
            "twentyeight",
            "clean17",
            "clean18",
            "clean19",
            "Low_15",
              "Low_16",
              "Low_17",
              "Low_18",
              "Low_19",
              "Low_20",
              "Low_21",
              "Low_22",
              "Low_23",
              "Low_24",
              "Low_25",
              "Low_26",
              "Low_27",
              "Low_28",
              "Dark_0",
              "Dark_3",
              "Dark_5",
              "Dark_6"]

    var_array = [zero,
            one,
            two,
            three,
            four,
            five,
            six,
            seven,
            eight,
            nine,
            ten,
            eleven,
            twelve,
            thirteen,
            fourteen,
            fifteen,
            sixteen,
            seventeen,
            eighteen,
            nineteen,
            twenty,
            twentyone,
            twentytwo,
            twentythree,
            twentyfour,
            twentyfive,
            twentysix,
            twentyseven,
            twentyeight,
            clean17,
            clean18,
            clean19,
             Low_15,
             Low_16,
             Low_17,
             Low_18,
             Low_19,
             Low_20,
             Low_21,
             Low_22,
             Low_23,
             Low_24,
             Low_25,
             Low_26,
             Low_27,
             Low_28,
             Dark_0,
             Dark_3,
             Dark_5,
             Dark_6]

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
              "red",
              "blue",
              "yellow",
              "green",
              "darkorange",
              "magenta",
              "cyan",
              "cyan",
              "red",
              "green",
              "green",
              "black",
              "yellow",
              "orange",
              "lime",
              "gray",
              "chocolate",
              "red",
              "blue",
              "yellow",
              "green",
              "darkorange",
              "magenta",
              "cyan",
              "green",
              "darkorange",
              "magenta",
              "cyan"]

    active = np.array([0, #good, possibly dark?
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #very bright
                       0, #good, dim
                       0, #good, dim
                       0, #really good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #good, dim
                       0, #START of clean data
                       0,
                       0,
                       1, #start of -166 temp data
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
                       0, # Start of Dark Rate Data
                       0,
                       0,
                       0],
                       dtype=bool)

    labels = np.array(labels)[active]
    var_array = np.array(var_array)[active]
    colors = np.array(colors)[active]
    nbins = 2500
    #fig, ax = plt.subplots(1, 1)
    #plot_SiPM_trace(ax, SBCcode.get_event(active_event, 0, "PMTtraces")["PMTtraces"])

    # plt.ioff()
    areas = []
    #separate array for clean areas
    clean_areas = []
    max_times = []
    n_runs =  len(var_array)
    plt.ioff()


    #arbitrary, not sure why these were set like this
    left_lim = 925 #regular values
    right_lim = 1150 #regular values
#    left_lim = 400 #dark rate values
#    right_lim = 1000 #dark rate values
    # left_lim = 400
    # right_lim = 2000
    for run_ix in range(n_runs):
        sub_areas = []
        sub_clean_areas = []
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


                #bring all values to 0
#                init_avg = np.average(yd_1[0:50]) #avg for dark rate
#                yd_1 -= init_avg

                #use yd_1
                #apply a low pass filter to data to eliminate noise
                order = 3
                fs = 1e9       # sample rate, Hz
                cutoff = 4e6 #4e6 desired cutoff frequency of the filter, Hz

                #loop through orders and plot with current order
                #for o in range(order):

                # # Get the filter coefficients so we can check its frequency response.
#                b, a = butter_lowpass(cutoff, fs, order)
                # #Filter the data, and plot both the original and filtered signals.
#                y = butter_lowpass_filter(yd_1, cutoff, fs, order)

                # Plot the frequency response.
                # w, h = freqz(b, a, worN=8000)
                # plt.subplot(2, 1, 1)
                # plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
                # plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
                # plt.axvline(cutoff, color='k')
                # plt.xlim(0, 0.5*fs)
                # plt.title("Lowpass Filter Frequency Response")
                # plt.xlabel('Frequency [Hz]')
                # plt.grid()


                # Demonstrate the use of the filter.
                # First make some data to be filtered.
                # T = 5.0         # seconds
                # n = int(T * fs) # total number of samples
                # t = np.linspace(0, T, n, endpoint=False)
                # # "Noisy" data.  We want to recover the 1.2 Hz signal from this.
                # data = np.sin(1.2*2*np.pi*t) + 1.5*np.cos(9*2*np.pi*t) + 0.5*np.sin(12.0*2*np.pi*t)

                #Plot the trace in the time domain
                #plt.subplot(3, 2, o+1)
                #use xd instead of t
                # plt.plot(xd, yd_1, 'b-', label='data')
                # #plt.plot(xd, y, 'g-', linewidth=2, label='filtered data')
                # plt.xlabel('Time [sec]')
                # #plt.title("Order" +str(o+1))
                # plt.grid()
                # plt.legend()
                # plt.show()

                #End order for-loop


                #plt.subplots_adjust(hspace=0.35)
                #plt.show()
                ###############




                #good_indices = (xd < right_lim) & (xd > left_lim)
                avg = np.average(yd_1[250:left_lim]) #regular avg
#                avg = np.average(y[250:right_lim]) #avg for dark rate
                #print("Avg:",avg)
                #pmt_area = scipy.integrate.trapz(yd_1[good_indices]-avg,
                #                                 dx=1e9*pmt_data["dt"][trig, 1])
                # if np.any(yd_0 > 0.15):
                #     pmt_area = -1.
                # else:
                pmt_area = np.sum(yd_1[left_lim:right_lim]-avg) #no filter applied
#                pmt_area = np.sum(y[left_lim:right_lim]-avg) #filter applied
#                pmt_area = np.sum(y[left_lim:right_lim]) #filter applied, no avg
                #pmt_clean_area = np.sum(y[left_lim:right_lim]-avg)
                #pmt_area = abs(np.sum(yd_0[left_lim:right_lim]-avg))

                #print(np.sum(yd_0))
                sub_areas.append(pmt_area)
                #sub_clean_areas.append(pmt_clean_area)
                sub_max_times.append(1e9 * pmt_data["dt"][trig, 0] * np.argmax(pmt_data["traces"][trig, 0, :]))

            # plotter = SiPMPlotter(pmt_data, runid=os.path.basename(active_event), ev=ev,
            #                       areas=sub_areas,
            #                       left=left_lim, right=right_lim)
            #
            # plotter.plot_SiPM_trace()
            # plt.show()
        areas.append(sub_areas)
        #clean_areas.append(sub_clean_areas)
        max_times.append(sub_max_times)



    #plt.hist(areas, bins=nbins, fill=False, color=colors, histtype="step", stacked=False, label=labels)
    (n, bins, patches) = plt.hist(areas, bins=nbins, fill=False, color=colors, histtype="step", stacked=False, label=labels)
#(n_clean, bins_clean, patches_clean) = plt.hist(clean_areas, bins=nbins, fill=False, color="blue", histtype="step", stacked=False, label=labels)


    #Attempt to Fit Gaussian Curve to each run
    if n_runs > 1:
        for num in range(n_runs): #aka number of sets of n
            run_data = n[num]

            #only use n and bins from desired domain
            #find indexes of bins from values lower_thresh to upper_thresh
            #        lower_thresh = .270 #.281
            #        upper_thresh = .684
            lower_thresh = np.inf
            peak = -1
            #        upper_thresh = -1
            start = 0
            peak_i = 0;
            end = 0

            #loop through the run data and find min and max values- range(nbins) is the same as looping through run_data
            #first, find the lower threshold
            #declare values to keep track of trend
            dec_stop_count = 0

            for i in range(nbins):
                d = run_data[i]
                if dec_stop_count > 50:
                    print("I",i)
                    break
                elif bins[i]>= 0 and d < lower_thresh and d > 0:
                    lower_thresh = d
                    start = i
                    dec_stop_count = 0
                elif bins[i]> 0:
                    dec_stop_count += 1

            #adjust start to the right an arbitrary amount
            #        start += 3
            #        lower_thresh = n[start]


            #then, find the peak
            for i in range(nbins-(start)):
                d = run_data[i+start]
                if d > peak and d > 0:
                    peak = d
                    peak_i = i+start

            #then calculate the upper threshold symmetric about the peak
            end = peak_i + (peak_i-start)
            #adjust end to the left by 25% of the difference
            end -= int(round((peak_i-start)/4))
            upper_thresh = run_data[end]

            print("\nGaussian Information for Dataset", labels[num])
            print("Start:",bins[start], "start value:", lower_thresh)
            print("Peak:",bins[peak_i], "Peak value:", peak)
            print("End:",bins[end], "end value:", upper_thresh)
            less_bins = bins[start:end]
            less_n = run_data[start:end]



            mean1 = 1
            sigma1 = .15
            popt,pcov = curve_fit(gaussian,less_bins,less_n,p0=[150,mean1,sigma1])

        #    plt.plot(less_bins,gaussian(less_bins,*popt),'ro:',label='fit')
            plt.plot(less_bins,gaussian(less_bins,*popt),'b-',label=labels[num]+' fit')
            #plot the lower thresh, peak, and upper thresh
            #plt.plot([bins[start], bins[peak_i], bins[end]],[lower_thresh,peak,upper_thresh], 'ro:')
            print("Amplitude:",popt[0],"Center:",popt[1],"Std Dev:",math.sqrt(popt[2]))

            #try to fit sum of Gaussians
            #find index of -0.2 in bins
            start_pos = np.nonzero(bins>-0.15)[0][0]
            #print("start pos", start_pos, "end",end)
            #get bin width
            bin_width = bins[1] - bins[0]
            #calculate upper thresh
            end = end-int(.01/bin_width)
            sum_bins = bins[start_pos:end] #from -0.2 to upper thresh - .12
            sum_n = run_data[start_pos:end]

            mean2 = 1
            sigma2 = .07
            popt_sum,pcov_sum = curve_fit(sum_gaussian,sum_bins,sum_n,p0=[150,mean1,sigma,1000,mean2,sigma2])
            plt.plot(sum_bins,sum_gaussian(sum_bins,*popt_sum),'g-',label=labels[num]+' sum fit')
            print("Noise Amplitude:",popt_sum[0],"Noise Center:",popt_sum[1],"Noise Std Dev",math.sqrt(popt_sum[2]/2),"1 Photon Amplitude:",popt_sum[3],"1 Photon Center:",popt_sum[4],"1 Photon Std Dev:",math.sqrt(popt_sum[5]/2))



    #else, there is only one run being run
    else:
        #only use n and bins from desired domain
        #find indexes of bins from values lower_thresh to upper_thresh
#        lower_thresh = .270 #.281
#        upper_thresh = .684
        lower_thresh = np.inf
        peak = -1
        #        upper_thresh = -1
        start = 0
        peak_i = 0;
        end = 0

        #loop through the run data and find min and max values- range(nbins) is the same as looping through run_data
        #first, find the lower threshold
        #declare values to keep track of trend
        dec_stop_count = 0

        for i in range(nbins):
            d = n[i]
            if dec_stop_count > 50:
                break
            elif bins[i]>= 0 and d < lower_thresh and d > 0:
                lower_thresh = d
                start = i
                dec_stop_count = 0
            elif bins[i]> 0:
                dec_stop_count += 1

        #adjust start to the right an arbitrary amount
#        start += 3
#        lower_thresh = n[start]


        #then, find the peak
        for i in range(nbins-(start)):
            d = n[i+start]
            if d > peak and d > 0:
                peak = d
                peak_i = i+start

        #then calculate the upper threshold symmetric about the peak
        end = peak_i + (peak_i-start)
        #adjust end to the left by 25% of the difference
        end -= int(round((peak_i-start)/4))
        upper_thresh = n[end]


#
#        for i in range(nbins):
#            if bins[i]>lower_thresh:
#                start = i
#                break
#        for i in range(nbins-(start+1)):
#            if bins[i+start] > upper_thresh:
#                end = i+start
#                break
        print("\nGaussian Information for Dataset", labels)
        print("start:",bins[start], "start value:", lower_thresh)
        print("Peak:",bins[peak_i], "Peak value:", peak)
        print("end:",bins[end], "end value:", upper_thresh)
        less_bins = bins[start:end]
        less_n = n[start:end]

        mean1 = 1
        sigma1 = .15
        # popt,pcov = curve_fit(gaussian,less_bins,less_n,p0=[150,mean1,sigma1])
        # #    plt.plot(less_bins,gaussian(less_bins,*popt),'ro:',label='fit')
        # plt.plot(less_bins,gaussian(less_bins,*popt),'b-',label=labels[0]+' fit')
        # #plot the lower thresh, peak, and upper thresh
        # plt.plot([bins[start], bins[peak_i], bins[end]],[lower_thresh,peak,upper_thresh], 'ro:')
        # print("Amplitude:",popt[0],"Center:",popt[1],"Std Dev:",math.sqrt(popt[2]))


        #try to fit sum of Gaussians
        #find index of -0.2 in bins
        start_pos = np.nonzero(bins>-0.15)[0][0]
        #print("start pos", start_pos, "end",end)
        #get bin width
        bin_width = bins[1] - bins[0]
        #calculate upper thresh
        end = end-int(.01/bin_width)
        sum_bins = bins[start_pos:end] #from -0.2 to upper thresh - .12
        sum_n = n[start_pos:end]

        mean2 = 1
        sigma2 = .07
        popt_sum,pcov_sum = curve_fit(sum_gaussian,sum_bins,sum_n,p0=[150,mean1,sigma1,1000,mean2,sigma2])
        plt.plot(sum_bins,sum_gaussian(sum_bins,*popt_sum),'r-',label=labels[0]+' Gaussian sum fit')
        print("Noise Amplitude:",popt_sum[0],"Noise Center:",popt_sum[1],"Noise Std Dev",math.sqrt(popt_sum[2]/2),"1 Photon Amplitude:",popt_sum[3],"1 Photon Center:",popt_sum[4],"1 Photon Std Dev:",math.sqrt(popt_sum[5]/2))

        #plot the individual gaussians
        #truncate popt_sum to noise and 1 photon fit
        popt_noise = popt_sum[0:3]
        popt_1 = popt_sum[3:6]
        #Calculate ranges for noise and single photon Gaussians
        noise_bins = np.arange(-.30,.30,bin_width)
        pixel_center = popt_sum[4]
        pixel_bins = np.arange(-.03,pixel_center+(pixel_center+.03),bin_width)
        plt.plot(noise_bins,gaussian(noise_bins,*popt_noise),'b-',label=labels[0]+' Noise fit')
        plt.plot(pixel_bins,gaussian(pixel_bins,*popt_1),'m-',label=labels[0]+' SiPM fit')

        #Calculate sigma for SiPM
        noise_var = popt_sum[2]/2
        SiPM_var = popt_sum[5]/2
        d_var = SiPM_var - noise_var
        SiPM_sigma = math.sqrt(d_var)
        print("SiPM Sigma", SiPM_sigma)


    plt.legend(loc="upper right")
    plt.suptitle("Histogram of Area")
    #plt.yscale("log")
    #plt.yscale("linear")
    plt.xlabel("Area")
    plt.ylabel("Count")

    #
    # plt.figure()
    # plt.hist(max_times, bins=nbins, fill=False, color=colors, histtype="step", stacked=False)
    # plt.legend(labels)
    # plt.suptitle("Time of max")
    # plt.ion()
    plt.grid()
    plt.show()




    pass
