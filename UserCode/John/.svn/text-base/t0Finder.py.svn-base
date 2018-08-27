## John Gresl
## 7/11/2018

import os
import numpy as np #
import scipy.signal
import matplotlib.pyplot as plt
import copy

from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE
from SBCcode.DataHandling.ReadBinary import ReadBlock as RB

import time

def extend_window(w, r):
    # Inputs:
    #   w: An array of 2 elements. Normally, this will be a window like [t1, t2]
    #   r: A float used as a ratio to extend w
    # Outputs: A rescaled version of w
    mp = 0.5*(w[1]+w[0])  ## Midpoint
    new_len = (w[1]-w[0])*(1+r) ## Length of new window
    return [mp-new_len/2, mp+new_len/2]

def freq_filter(freqs, lower=None, upper=None):
    # Inputs:
    #   freqs: An array of frequency bins
    #   lower: The lower frequency to cut-off at
    #   upper: The upper frequency to cut-off at
    # Outputs: An array of indices where the the frequency in freqs is between lower and upper
    if lower is None and upper is None:
        return freqs
    if lower is None:
        return np.where([x <=  upper for x in freqs])
    if upper is None:
        return np.where([x >= lower for x in freqs])
    return np.where([lower <= x <= upper for x in freqs])


def spectrum_sums2(spectrum, fr, N, lowerf=None, upperf=None):
    # Inputs:
    #   spectrum: The output 2d spectrum from a spectogram
    #   fr: A list of frequency bins corresponding to the spectrum
    #   N: Number of bins
    #   lowerf: The lower frequency to cut-off at
    #   upperf: The upper frequency to cut-off at
    # Outputs: A compressed 1d array where each element is the sum of a bin from spectrum, only counting
    #          frequencies between lowerf and upperf
    out = []
    good_indices = freq_filter(fr, lowerf, upperf)
    for n in range(N):
        out.append(np.trapz(spectrum[good_indices[0], n], dx = np.mean(np.diff(fr))))
        #out.append(sum(spectrum[good_indices[0], n]))
    return out


def rescale_window(w1, w2):
    # Inputs:
    #   w1: An array with 2 elements
    #   w2: An array with 2 elements
    # Outputs: A rescaled version of w2 so tha the endpoints of w2 match w1 but the number of elements remain the same
    y1, y2 = min(w1), max(w1)
    x1, x2 = min(w2), max(w2)
    if x1==x2:
        return 0*w2
    a = (y1-y2)/(x1-x2)
    b = (x1*y2-x2*y1)/(x1-x2)
    return a*w2+b


def corr_signal(tau, dt, t0, N, fit_type = 0, shift = 10):
    # Inputs:
    #   tau: Time constant on exponential decay
    #   dt: Step size for the x-axis
    #   t0: Where the exponential signal will start. Not important when used with correlation
    #   N: Number of points requested
    #   fit_type: The type of signal to create. See corr_signal_type_templates.txt for a better explanation.
    #               fit_type = 0 --> Exponential decay
    #               fit_type = 1 --> Constant 1 followed by exponential decay (continuous)
    #               fit_type = 2 --> Linear increase followed by exponential decay
    #               fit_type = 3 --> Log increase followed by exponential decay
    #               fit_type = 4 --> 0 value followed by an exponential decrease. Discontinuous.
    # Outputs:
    #   t: t-values for plotting
    #   y: y-values of our filter signal.
    # After careful analysis, we've determined that there reaches a point in the filtered piezo signal that
    # exhibits a sharp increase followed by an exponential decay. This function returns a brief exponential
    # decay function for use with convolution/correlation.
    shift = int(np.ceil(shift))
    t = np.linspace(t0, t0+dt*N, N)
    y = np.exp(-(t-t0)/tau)
    ycopy = copy.deepcopy(y)
    if fit_type == 0:
        pass
    elif fit_type == 1:
        for n in range(len(y) - shift):
            y[n+shift] = ycopy[n]
        y[0:shift] = 1
    elif fit_type == 2:
        for n in range(len(y) - shift):
            y[n + shift] = ycopy[n]
        y[0:shift] = (t[0:shift] - t0)/(shift*dt)
    elif fit_type == 3:
        for n in range(len(y) - shift):
            y[n + shift] = ycopy[n]
            y[0:shift] = np.log((t[0:shift] + 1 - t0)) / np.log(shift*dt + 1)
    elif fit_type == 4:
        for n in range(len(y) - shift):
            y[n+shift] = ycopy[n]
            y[0:shift] = 0
    return t, y

#t,y = corr_signal(1.5, 0.1, 0, 45, 3, 20) # <-- Uncomment to test different fit types
#plt.ioff()
#plt.plot(t,y)
#plt.show()

def find_t0_from_corr(corrt, corry):
    # Inputs:
    #   corrt: Time-values of the correlation signal
    #   corry: Y-values of the correlation signal
    # Outputs: The time of the maximum in corry such that corrt is less than or equal to 0.
    n = np.where(corrt>=0)
    corry[n]=0
    return corrt[np.argmax(corry)]


def run(data, tau, view_plots=False):
    LOWER = 20000
    UPPER = 40000
    try:
        piezo1 = data["fastDAQ"]["Piezo1"]
        piezo2 = data["fastDAQ"]["Piezo2"]
        timebase = data["fastDAQ"]["time"]
        xextent = (min(timebase), max(timebase))
        dt = np.mean(np.diff(timebase))
    except:
        print("Unable to load piezo data. Returning.")
        return
    if view_plots:
        fig, ax = plt.subplots(nrows=5, ncols=2, sharex=True)

        spec1_ax = ax[0][0] ; spec2_ax = ax[0][1]
        spec1_ax.set_title("Piezo1") ; spec2_ax.set_title("Piezo2")
        raw1_ax = ax[1][0] ; raw2_ax = ax[1][1]
        filt1_ax = ax[2][0] ; filt2_ax = ax[2][1]
        corr1_ax = ax[3][0] ; corr2_ax = ax[3][1]
        stat1_ax = ax[-1][0] ; stat2_ax = ax[-1][1]

        sp1, fr1, bn1, im1 = spec1_ax.specgram(piezo1, Fs=1./dt, NFFT=512, noverlap=450, xextent=xextent, mode="psd")
        sp2, fr2, bn2, im2 = spec2_ax.specgram(piezo2, Fs=1./dt, NFFT=512, noverlap=450, xextent=xextent, mode="psd")

        spec1_ax.set_ylim(extend_window([LOWER, UPPER], 0.05))
        spec2_ax.set_ylim(extend_window([LOWER, UPPER], 0.05))

        raw1_ax.plot(timebase, piezo1)
        raw2_ax.plot(timebase, piezo2)

        N1 = len(bn1)
        N2 = len(bn2)

        sp1_sums = spectrum_sums2(sp1, fr1, N1, LOWER, UPPER)
        sp2_sums = spectrum_sums2(sp2, fr2, N2, LOWER, UPPER)
        sp1_sums = scipy.signal.medfilt(sp1_sums)
        sp2_sums = scipy.signal.medfilt(sp2_sums)


        rescaled_t1 = rescale_window(xextent, bn1)
        rescaled_t2 = rescale_window(xextent, bn2)

        filt1_ax.plot(rescaled_t1, sp1_sums)
        filt2_ax.plot(rescaled_t2, sp2_sums)

        corr_dt = np.mean(np.diff(rescaled_t1))
        corr_n = 1000

        #corr_signal(tau, dt, t0, N, fit_type, shift)
        corr1_signalx, corr1_signaly = corr_signal(tau, corr_dt, rescaled_t1[0], corr_n, fit_type=2, shift=corr_n/100)
        corr2_signalx, corr2_signaly = corr_signal(tau, corr_dt, rescaled_t2[0], corr_n, fit_type=2, shift=corr_n/100)

        corr1 = np.correlate(sp1_sums, corr1_signaly, "same")
        corr2 = np.correlate(sp2_sums, corr2_signaly, "same")

        corr1t = rescaled_t1 - 0.5 * corr_n * corr_dt
        corr2t = rescaled_t2 - 0.5 * corr_n * corr_dt
        corr1_ax.plot(corr1t, corr1, "b-")

        corr2_ax.plot(corr2t, corr2, "b-")


        # reco_data = RB(os.path.join(RECO_DIRECTORY, "{}/AcousticAnalysis_{}.bin".format(event_id, event_id)))
        # old_t0 = sum(reco_data["bubble_t0"][int(event)]) / 2

        new_t01 = find_t0_from_corr(corr1t, corr1)
        new_t02 = find_t0_from_corr(corr2t, corr2)

        sub0ind1 = np.where(corr1t < 0)
        sub0ind2 = np.where(corr2t < 0)
        sub0sp1_sums = sp1_sums[sub0ind1]
        sub0sp2_sums = sp2_sums[sub0ind2]

        filt1_ax.plot(corr1_signalx + (new_t01-rescaled_t1[0]), 0.5*max(sub0sp1_sums) * corr1_signaly, "g-", linewidth=2)
        filt2_ax.plot(corr2_signalx + (new_t02-rescaled_t2[0]), 0.5*max(sub0sp2_sums) * corr2_signaly, "g-", linewidth=2)
        n = len(ax.flatten())
        for subaxn in range(n):
            if subaxn%2 == 0:
                ax.flatten()[subaxn].axvline(x=new_t01, color="g", linewidth=1.5)
                pass
            else:
                ax.flatten()[subaxn].axvline(x=new_t02, color="g", linewidth=1.5)
            ax.flatten()[subaxn].set_xlim(xextent)
        plt.show()

    if not view_plots: # Then don't waste time generating them
        t0 = time.time()
        fr1, bn1, sp1 = scipy.signal.spectrogram(piezo1, fs=1./dt, nfft=512, noverlap=450, mode="psd", window="hanning", nperseg=512)
        fr2, bn2, sp2 = scipy.signal.spectrogram(piezo2, fs=1./dt, nfft=512, noverlap=450, mode="psd", window="hanning", nperseg=512)
        t1 = time.time()
        print("\tSpecgram took {}s".format(t1-t0))
        N1 = len(bn1)
        N2 = len(bn2)

        sp1_sums = spectrum_sums2(sp1, fr1, N1, LOWER, UPPER)
        sp2_sums = spectrum_sums2(sp2, fr2, N2, LOWER, UPPER)
        sp1_sums = scipy.signal.medfilt(sp1_sums)
        sp2_sums = scipy.signal.medfilt(sp2_sums)

        rescaled_t1 = rescale_window(xextent, bn1)
        rescaled_t2 = rescale_window(xextent, bn2)

        corr_dt1 = np.mean(np.diff(rescaled_t1))
        corr_dt2 = np.mean(np.diff(rescaled_t2))
        corr_n1 = 1000
        corr_n2 = 1000

        corr1_signalt, corr1_signaly = corr_signal(tau, corr_dt1, rescaled_t1[0], corr_n1, fit_type=0)
        corr2_signalt, corr2_signaly = corr_signal(tau, corr_dt2, rescaled_t2[0], corr_n2, fit_type=4, shift=corr_n2/150)

        corr1 = np.correlate(sp1_sums, corr1_signaly, "same")
        corr2 = np.correlate(sp2_sums, corr2_signaly, "same")
        t2 = time.time()
        print("\tCorr took {}s".format(t2-t1))
        corr1t = rescaled_t1 - 0.5 * corr_n1 * corr_dt1
        corr2t = rescaled_t2 - 0.5 * corr_n2 * corr_dt2

        new_t01 = find_t0_from_corr(corr1t, corr1)
        new_t02 = find_t0_from_corr(corr2t, corr2)


    return new_t01, new_t02


if __name__ == "__main__":
    TAU_PEAK = 0.0025884277467056165  # <-- This comes from TauResultAnalysis.py
    TAU_AVERAGE = 0.0038163479219674467 # <-- This also ^^
    RAW_DIRECTORY = "/bluearc/storage/SBC-17-data/"
    RECO_DIRECTORY = "/pnfs/coupp/persistent/grid_output/SBC-17/output"
    VIEW_PLOTS = True
    plt.ioff()
    for i in range(1):
        data = GE(os.path.join(RAW_DIRECTORY, "20170928_5"), 36, "fastDAQ")
        t0 = time.time()
        print(run(data, TAU_AVERAGE, VIEW_PLOTS))
        t1 = time.time()
        print("Event: {}. Finding t0 took {} seconds".format(i, t1-t0))

    #run(RAW_DIRECTORY, "20170630_0", "0", TAU_AVERAGE, VIEW_PLOTS)
    #run(RAW_DIRECTORY, "20170627_0", "3", TAU_AVERAGE, VIEW_PLOTS)
