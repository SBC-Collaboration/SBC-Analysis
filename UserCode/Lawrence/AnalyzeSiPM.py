import sys
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

def process_input(argv):
    #User Input Form: "python AnalyzeSiPM.py [run_file_path] [hist_bool] [filter_order] [traces_bool] [channel_num]"
    #First, parse filepath
    try:
        filepath = sys.argv[1]
        print ("Analyzing SiPM Data at filepath:", filepath)
    except:
        sys.exit("Error: No filepath provided. Accepted format is:\n python AnalyzeSiPM.py [run_file_path] [graph_bool=1] [filter_order=0] [traces_bool=0] [channel_num=0]")
    #Start Processing Options
    print ("\nOptions:")
    #Next, get hist_bool
    try:
        hist_bool = int(sys.argv[2])
        if hist_bool == 0:
            print ("Plot Histogram OFF") #For debugging
        else: #if hist_bool is anything except 0, hist will be turned on
            hist_bool = 1
            print ("Plot Histogram ON") #For debugging
    except:
        #If no input, the program should plot the histogram
        hist_bool = 1
        print ("Plot Histogram ON") #For debugging

    #Next, get filter_order
    try:
        filter_order = int(sys.argv[3])
        if filter_order <= 5 and filter_order > 0:
            print ("Filter ON with order", filter_order) #For debugging
        else: #if filter_order is anything out of range 1-5 inclusive, filter will be turned off
            filter_order = 0
            print ("Filter OFF") #For debugging
    except:
        #If no input, don't filter
        filter_order = 0
        print ("Filter OFF") #For debugging

    #Lastly, get traces_bool
    try:
        traces_bool = int(sys.argv[4])
        if traces_bool == 1:
            print ("Graph Traces ON") #For debugging
        else: #if traces_bool is anything except 1, traces graphing will be turned off
            traces_bool = 0
            print ("Graph Traces OFF") #For debugging
    except:
        #If no input, don't plot traces
        traces_bool = 0
        print ("Graph Traces OFF") #For debugging

    #Next, get the channel number
    try:
        channel_num = int(sys.argv[5])
        if channel_num == 0:
            print ("Processing Data from Channel 0") #For debugging
        else: #if channel_num is anything except 0, channel 1 will be processed
            channel_num = 1
            print ("Processing Data from Channel 1") #For debugging
    except:
        #If no input, process channel 0
        channel_num = 1
        print ("Processing Data from Channel 1") #For debugging

    return [filepath, hist_bool, filter_order, traces_bool, channel_num]

def FilterData(y, filter_order):
    #Arrays to keep track of data
    sub_clean_areas = []

    #Low Pass Filter Params
    fs = 1e9        #Sample rate (Hz)
    cutoff = 4e6    #Filter cutoff frequency (Hz)

    # Get the filter coefficients so we can check its frequency response if we want
    b, a = butter_lowpass(cutoff, fs, filter_order)
    #Filter the data and return
    cleaned_y = butter_lowpass_filter(y, cutoff, fs, filter_order)
    return cleaned_y

def GetAreas(filepath, channel_num, filter_order):
    #Parameters for calculating area of the traces
    left_lim = 925      #Left limit index value of timestep to begin calculating area
    right_lim = 1150    #Right limit index value of timestep to begin calculating area
    avg_start = 250     #The index value of timestep to begin calculation of average

    #arrays to hold data
    areas = []
    max_times = []
    time_arr = []
    y_arr = []
    clean_y_arr = []

    #Get events from run file
    events = SBCtools.BuildEventList(filepath)

    #Loop through events
    for ev in events:
        #Get traces from event
        pmt_data = SBCcode.get_event(filepath, ev, "PMTtraces", max_file_size=1300)["PMTtraces"]
        n_triggers = pmt_data["traces"].shape[0]

        #Loop through traces
        for trig in range(n_triggers):
            #Calculate the time value
            xd = 1e9 * np.arange(pmt_data["traces"].shape[2]) * pmt_data["dt"][trig, 0]
            #Calculate voltage value from the correct channel
            if channel_num == 0:
                y = pmt_data["traces"][trig, 0, :] * pmt_data["v_scale"][trig, 0] + \
                       pmt_data["v_offset"][trig, 0]
            else:
                y = pmt_data["traces"][trig, 1, :] * pmt_data["v_scale"][trig, 1] + \
                   pmt_data["v_offset"][trig, 1]

            #Filter the noise if the order is valid
            if filter_order != 0:
                cleaned_y = FilterData(y, filter_order)
                avg = np.average(cleaned_y[avg_start:left_lim]) #Cleaned Avg
                pmt_area = np.sum(cleaned_y[left_lim:right_lim]-avg) #Calculate Area
                #store clean y
                clean_y_arr.append(cleaned_y)

            #If no filter applied
            else:
                avg = np.average(y[avg_start:left_lim]) #Average of raw data
                pmt_area = np.sum(y[left_lim:right_lim]-avg) #Calculate Area

            #Append area to areas array
            areas.append(pmt_area)

            #Store time and y arrays
            y_arr.append(y)
            time_arr.append(xd)

            # #Append max times to array. Not sure if necessary, but will keep here for now
            # max_times.append(1e9 * pmt_data["dt"][trig, 0] * np.argmax(pmt_data["traces"][trig, 0, :]))

    #Return values
    return [areas, time_arr, y_arr, clean_y_arr]

def FindSinglePeakBounds(n, nbins):
    #lower_thresh should be xval ~.270
    #upper_thresh should be xval ~.684

    #Variable values to be obtained
    lower_thresh = np.inf
    upper_thresh = 0
    peak = -1
    start = 0
    peak_i = 0
    end = 0
    #Decrement Stop Count to keep track of when the data stops decreasing
    dec_stop_count = 0

    #Loop through run data
    #First, find the lower threshold
    for i in range(nbins):
        #Get current histogram val
        data = n[i]
        #If the hist data has not decreased for 50 timesteps, stop looping
        if dec_stop_count > 50:
            break
        #If new min found, adjust lower_thresh. Only care about positive x vals
        elif bins[i]>= 0 and data < lower_thresh and data > 0:
            lower_thresh = data
            start = i
            dec_stop_count = 0
        elif bins[i]> 0:
            dec_stop_count += 1

    #Next, find the peak
    #This loops through all the remainig data, so it's kind of inefficient
    #Could use an increment_stop_count so that it stops looping once 50 timesteps pass without finding an increase
    for i in range(nbins-(start)):
        d = n[i+start]
        if d > peak and d > 0:
            peak = d
            peak_i = i+start

    #Next, calculate the upper threshold symmetric about the peak
    end = peak_i + (peak_i-start)

    #adjust end to the left by 25% of the difference
    end -= int(round((peak_i-start)/4))
    upper_thresh = n[end]

    return [start,peak_i,end,lower_thresh,peak,upper_thresh]

def FitGauss(n, bins, nbins):
    #only use n and bins from desired domain
    #Process Single-Pixel peak: find index of end xval
    [start,peak_i,end,lower_thresh,peak,upper_thresh] = FindSinglePeakBounds(n, nbins)

    #Initial guess parameters for curve fit
    amp_noise = 1000
    mean_noise = 1
    sigma_noise = .15
    amp_SP = 150
    mean_SP = 1
    sigma_SP = .07

    #Fit Sum of Gaussiants
    #First, find index of -0.15 in bins. This
    start_pos = np.nonzero(bins>-0.15)[0][0]
    #Get bin width
    bin_width = bins[1] - bins[0]
    #Re-adjust upper thresh-- Without this, the curve fit fails
    end = end-int(.01/bin_width)
    #Truncate data to only hold relevant data in the desired domain
    sum_bins = bins[start_pos:end] #from -0.15 to upper thresh
    sum_n = n[start_pos:end]

    #Fit the Sum of Gaussians
    popt_sum,pcov_sum = curve_fit(sum_gaussian,sum_bins,sum_n,p0=[150,mean_noise,sigma_noise,1000,mean_SP,sigma_SP])

    # #This is for visualizing the range of data being approximated to the Gaussian Curve (avg_start, left_lim, right_lim)
    # #Uncomment this to plot these points on the histogram
    #ax.plot([bins[start_pos], bins[start], bins[peak_i], bins[end]],[n[start_pos],lower_thresh,peak,upper_thresh], 'ro')


    #Return the Gaussian Params
    return [popt_sum, sum_n, sum_bins, bin_width, start]

def PlotIndividualGaussians(popt_sum, bin_width):
    #Truncate popt_sum to Noise and Single Photon fit
    popt_noise = popt_sum[0:3]
    popt_SP = popt_sum[3:6]

    #Calculate ranges for noise and single photon Gaussians
    noise_bins = np.arange(-.30,.30,bin_width)
    pixel_center = popt_sum[4]
    pixel_bins = np.arange(-.03,pixel_center+(pixel_center+.03),bin_width)
    ax.plot(noise_bins,gaussian(noise_bins,*popt_noise),'b--',label='Noise fit')
    ax.plot(pixel_bins,gaussian(pixel_bins,*popt_SP),'g--',label='SiPM fit')
    return

def CalculateSiPMSigma(popt_sum):
    #Calculate variance of noise
    noise_var = popt_sum[2]/2
    #Calculate variance of SP peak
    SiPM_var = popt_sum[5]/2
    #Calculate Difference
    d_var = abs(SiPM_var - noise_var)
    #Convert variance to std. dev.
    SiPM_sigma = math.sqrt(d_var)
    return SiPM_sigma


def PlotTraces(n,bins,local_min,nbins,bin_width,time_arr, y_arr, clean_y_arr):
    print ("\nPlotting Traces")
    #First, find max of n
    noise_cen = bins[np.argmax(n)]
    SP_cen = bins[np.argmax(n[local_min:])] #this could be the problem

    #Find first 10 traces that contribute to noise peak
    #and first 10 traces that contribute to single pixel peak
    noise_indexes = []
    SP_indexes = []
    noise_count = 0
    SP_count = 0
    for i in range(nbins): #loop through areas array
        if areas[i] > SP_cen and areas[i] < SP_cen + bin_width and SP_count < 10: #find first 10 instances of SP within bin
            SP_indexes.append(i)
            SP_count += 1
        elif areas[i] > noise_cen and areas[i] < noise_cen + bin_width and noise_count < 10: #find first 10 instances of noise within bin
            noise_indexes.append(i)
            noise_count += 1
        elif noise_count >= 10 and SP_count >= 10:
            break

    print ("Noise indexes", noise_indexes)
    print ("SP indexes", SP_indexes)


    #Get traces and store
    noise_traces = []
    SP_traces = []
    #Loop through 10 indexes
    for i in range(10):
        curr_noise_i = noise_indexes[i]
        curr_SP_i = SP_indexes[i]
        noise_traces.append(y_arr[curr_noise_i])
        SP_traces.append(y_arr[curr_SP_i])

    #Plot traces in the time domain
    plt.figure()
    #Plot Noise
    plt.subplot(2,1,1)
    plt.plot(time_arr, noise_traces, 'b-', label='Noise')
    #Plot filtered noise
    if (filter_order != 0):
        cleaned_noise = []
        for i in range(10):
            curr_noise_i = noise_indexes[i]
            cleaned_noise.append(clean_y_arr[curr_noise_i])
        plt.plot(time_arr, cleaned_noise, 'g-', linewidth=2, label='Filtered Noise')
    #Plot SP Pulses
    plt.subplot(2,1,2)
    plt.plot(time_arr, y, 'b-', label='Single Photon Pulse')
    #Plot filtered SP pulses
    if (filter_order != 0):
        cleaned_pulses = []
        for i in range(10):
            curr_SP_i = SP_indexes[i]
            cleaned_pulses.append(clean_y_arr[curr_SP_i])
        plt.plot(time_arr, cleaned_pulses, 'g-', linewidth=2, label='Filtered Single Photon Pulse')
    plt.xlabel('Time (sec)')
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    #Parse arguments from user input
    [filepath, hist_bool, filter_order, traces_bool,channel_num] = process_input(sys.argv)

    #Non-User Defined Parameters
    nbins = 2500

    #get trace areas
    [areas, time_arr, y_arr, clean_y_arr] = GetAreas(filepath, channel_num, filter_order)

    #Plot the data on a histogram
    fig,ax = plt.subplots(1)
    (n, bins, patches) = ax.hist(areas, bins=nbins, fill=False, color=["magenta"], histtype="step", stacked=False, label=filepath)

    #Fit Sum of Gaussians and get params
    [popt_sum, sum_n, sum_bins,bin_width,local_min_i] = FitGauss(n, bins, nbins)

    #Print the Parameters
    print("\n#DATA#\nNoise Amplitude:")
    print(popt_sum[0],"\nNoise Center:")
    print(popt_sum[1],"\nNoise Std Dev:")
    print(math.sqrt(popt_sum[2]/2),"\nSingle Photon Amplitude:")
    print(popt_sum[3],"\nSingle Photon Center:")
    print(popt_sum[4],"\nSingle Photon Std Dev:")
    print(math.sqrt(popt_sum[5]/2))

    #Calculate sigma for SiPM
    SiPM_sigma = CalculateSiPMSigma(popt_sum)
    print("SiPM Sigma:")
    print(SiPM_sigma)

    #Show Histogram
    if hist_bool:
        #Plot the Gaussian Sum on the histogram
        ax.plot(sum_bins,sum_gaussian(sum_bins,*popt_sum),'r-',label='Gaussian Sum Fit')

        #Plot the individual Gaussians
        PlotIndividualGaussians(popt_sum,bin_width)

        #Plot the information on the figure
        fig.subplots_adjust(bottom=0.4)
        textstr = "Noise Amplitude: " + str(popt_sum[0]) + "\nNoise Center: " + str(popt_sum[1]) + "\nNoise Std Dev: " + str(math.sqrt(popt_sum[2]/2)) + "\nSingle Photon Amplitude: " + str(popt_sum[3]) + "\nSingle Photon Center: " + str(popt_sum[4]) + "\nSingle Photon Std Dev:" + str(math.sqrt(popt_sum[5]/2)) + "\n SiPM Sigma:" + str(SiPM_sigma)
        fig.text(.5, .05, textstr, ha='center')

        #Finish Plotting
        plt.legend(loc="upper right")
        plt.suptitle("Histogram of Pulse Area")
        plt.yscale("log")
        plt.xlabel("Area")
        plt.ylabel("Count")
        plt.grid()
        plt.show()

    #Plot traces
    # if traces_bool:
    #     PlotTraces(n,bins,local_min_i,nbins,bin_width, time_arr, y_arr, clean_y_arr)
    pass
