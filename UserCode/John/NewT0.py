'''
New t0 algorithm is getting close, but not close enough.

1. Run previous analysis using only exponential signal. In general, this should predict a t0 that is slightly ahead of the actual t0.
2. Determine a baseline by looking inbetween the LED pulses and taking the average. (Integrate and divide by length?)
3. From the previous t0, start looking backwards in the filtered signal (the signal before correlations) until we find a region that is below the baseline for at least 100 ms.
	a. The new t0 should be the right end-point of this region plus one half the width of our FFT window from the spectrogram

'''

import numpy as np
import scipy.signal
import copy
import operator
import re
# import matplotlib.pyplot as plt
import copy
# from SBCcode.Tools import SBCtools
from ...DataHandling.WriteBinary import WriteBinaryNtupleFile as WB

def get_runs(dir, search_for="folders", do_sort=True, reverse=False):
    # Author: John Gresl
    # Inputs:
    #   dir: A directory containing folders of all runs
    #   search_for: A string indicating whether to extract run_ids from folders or from files
    #               Can be either "folders" or "files"
    #   do_sort: A boolean. If true, will return the sorted run_list from sort_runs above.
    #   reverse: Bool. If true, returns the reverse-sorted list.
    # Outputs: An array of run-strings from a given directory. This is passed to sort_runs if do_sort is true
    if search_for.lower().strip() == "folders":
        base_dirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    elif search_for.lower().strip() == "files":
        base_dirs = [d for d in os.listdir(dir) if os.path.isfile(os.path.join(dir, d))]
    else:
        raise ValueError("'search_for' must be either 'folders' or 'files. Got {}".format(search_for))
    to_match = re.compile("([0-9]{8}_[0-9]+){1}")
    out = []
    for d in base_dirs:
        match_q = re.search(to_match, d)
        if match_q:
            out.append(match_q.group(0))
    if do_sort:
        return sort_runs(out, reverse=reverse)
    return out


def sort_runs(arr, reverse=False):
    # Author: John Gresl 8/2/2018
    # Input:
    #   arr: An array of run_ids as strings. Should look like ["20170623_0", "20170623_5", etc...]
    #   reverse: Bool. If true, returns the reverse-sorted list.
    # Outputs: A natural-ish sorted version that puts the dates in order and the run numbers for each date in order
    s = sorted([np.int32(runid.split("_")) for runid in arr], key=operator.itemgetter(0, 1), reverse=reverse)
    return ["_".join(np.str(runid).strip("[]").split()) for runid in s]


def trim_runlist(arr, start=None, stop=None):
    # Author: John Gresl 8/2/2018
    # Inputs:
    #   arr: An array of run_ids as strings. Should look like ["20170623_0", "20170623_5", etc...]
    #   start: Start run number. If this is not supplied, will start at the beginning
    #   stop: Stop run number. If this is not supplied, will continue to end
    # Outputs: A sorted, trimmed runlist that goes from start to stop
    arr = sort_runs(arr)
    start = arr[0] if start == None else start
    stop = arr[-1] if stop == None else stop
    start_date = int(start.split("_")[0])
    start_run_num = int(start.split("_")[1])
    stop_date = int(stop.split("_")[0])
    stop_run_num = int(stop.split("_")[1])
    out = []
    for run in arr:
        date = int(run.split("_")[0])
        run_num = int(run.split("_")[1])
        if start_date > date or date > stop_date:
            continue
        if (start_run_num > run_num and date == start_date) or (run_num > stop_run_num and date == stop_date):
            continue
        out.append(run)
    return out


def extend_window(w, r):
    # Inputs:
    #   w: An array of 2 elements. Normally, this will be a window like [t1, t2]
    #   r: A float used as a ratio to extend w
    # Outputs: A rescaled version of w
    mp = 0.5*(w[1]+w[0])  # Midpoint
    new_len = (w[1]-w[0])*(1+r)  # Length of new window
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
        return np.where([x <= upper for x in freqs])
    if upper is None:
        return np.where([x >= lower for x in freqs])
    return np.where([lower <= x <= upper for x in freqs])


def closest_index(arr, el):
    # Inputs:
    #   arr: A 1-dimensional array
    #   el: Any element
    # Outputs: The FIRST index of the item in arr that is closest to el.
    # Notes: Arr does NOT have to be sorted.
    return np.argmin(np.abs(arr-el))


def spectrum_sums(spectrum, fr, n, lowerf=None, upperf=None):
    # Inputs:
    #   spectrum: The output 2d spectrum from a spectogram
    #   fr: A list of frequency bins corresponding to the spectrum
    #   n: Number of bins
    #   lowerf: The lower frequency to cut-off at
    #   upperf: The upper frequency to cut-off at
    # Outputs: A compressed 1d array where each element is the sum of a bin from spectrum, only counting
    #          frequencies between lowerf and upperf
    out = []
    good_indices = freq_filter(fr, lowerf, upperf)
    for subn in range(n):
        out.append(np.trapz(spectrum[good_indices[0], subn], dx=np.mean(np.diff(fr))))
    return out


def rescale_window(w1, w2):
    # Inputs:
    #   w1: An array with 2 elements
    #   w2: An array with 2 elements
    # Outputs: A rescaled version of w2 so tha the endpoints of w2 match w1 but the number of elements remain the same
    y1, y2 = min(w1), max(w1)
    x1, x2 = min(w2), max(w2)
    if x1 == x2:
        return 0*w2
    a = (y1-y2)/(x1-x2)
    b = (x1*y2-x2*y1)/(x1-x2)
    return a*w2+b


def corr_signal(tau, dt, t0, n, fit_type=0, shift=10):
    # Inputs:
    #   tau: Time constant on exponential decay
    #   dt: Step size for the x-axis
    #   t0: Where the exponential signal will start. Not important when used with correlation
    #   N: Number of points requested
    #   fit_type: The type of signal to create. See corr_signal_type_templates.py for a better explanation.
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
    t = np.linspace(t0, t0+dt*n, n)
    y = np.exp(-(t-t0)/tau)
    ycopy = copy.deepcopy(y)
    if fit_type == 0:
        pass
    elif fit_type == 1:
        for subn in range(len(y) - shift):
            y[subn+shift] = ycopy[subn]
        y[0:shift] = 1
    elif fit_type == 2:
        for subn in range(len(y) - shift):
            y[subn + shift] = ycopy[subn]
        y[0:shift] = (t[0:shift] - t0)/(shift*dt)
    elif fit_type == 3:
        for subn in range(len(y) - shift):
            y[subn + shift] = ycopy[subn]
            y[0:shift] = np.log((t[0:shift] + 1 - t0)) / np.log(shift*dt + 1)
    elif fit_type == 4:
        for subn in range(len(y) - shift):
            y[subn+shift] = ycopy[subn]
            y[0:shift] = 0
    return t, y



def find_t0_from_corr(corrt, corry):
    # Inputs:
    #   corrt: Time-values of the correlation signal
    #   corry: Y-values of the correlation signal
    # Outputs: The time of the maximum in corry such that corrt is less than or equal to 0.
    n = np.where(corrt >= 0)
    corry[n] = 0
    return corrt[np.argmax(corry)]

def within_baseline(arr, baseline, rms):
    return np.all(arr<(baseline + 5*rms))

def calculate_t0(data, tau, lower, upper, piezo1_fit_type=0, piezo2_fit_type=4, view_plots=False):
    # Inputs:
    #   data: Data returned from SBCcode's GetEvent function. Must have fastDAQ loaded.
    #         GetEvent is found within SBCcode/DataHandling/GetSBCEvent.py
    #   tau: The time constant we are trying to fit to the exponential decay that occurs
    #        immediately after the bubble forms
    #   lower: The lower frequency threshold for cutting off the spectrogram
    #   upper: The upper frequency threshold for cutting off the spectrogram
    #   piezo1_fit_type: The type of fit to use when trying to match the filtered piezo1 signal. Defaults to 0.
    #   piezo2_fit_type: The type of fit to use when trying to match the filtered piezo2 signal. Defaults to 4.
    #                    For a description of fit_types, see corr_signal above or check corr_signal_types.py
    #   view_plots: Boolean. If true, will display some plots for analysis.
    # Outputs: A dictionary of results for the Acoustic Analysis.
    # Issues:
    #   1. The actual run_id and ev cannot be extracted from the GetEvent output, which would mean in addition
    #      to previously loading the GetEvent result, we also have to pass the event_id and ev to run() which
    #      in itself is fine, but it's bad programming practice. If there's a random GetEvent() output sitting
    #      around, we should be able to extract the event id and event from it...

    default_output = {"bubble_t0": np.zeros(2)-1.}  # Default output incase the file can't be opened or something
    try:
        if not data["fastDAQ"]["loaded"]:
            return default_output
        out = default_output
        piezo1 = data["fastDAQ"]["Piezo1"]
        piezo2 = data["fastDAQ"]["Piezo2"]
        #
        #
        # filter_d = np.exp(-.1)
        # piezo1 = scipy.signal.lfilter(np.float64([1])-filter_d, np.float64([1, -filter_d]), piezo1)
        # piezo2 = scipy.signal.lfilter(np.float64([1]) - filter_d, np.float64([1, -filter_d]), piezo2)
        timebase = data["fastDAQ"]["time"]
        xextent = (min(timebase), max(timebase))
        dt = np.mean(np.diff(timebase))
    except Exception as e:
        print("Unable to load piezo data.")
        return default_output

    try:
        if view_plots:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(nrows=4, ncols=2, sharex=True)
            spec1_ax = ax[0][0] ; spec2_ax = ax[0][1]
            spec1_ax.set_title("Piezo1") ; spec2_ax.set_title("Piezo2")
            raw1_ax = ax[1][0] ; raw2_ax = ax[1][1]
            filt1_ax = ax[2][0] ; filt2_ax = ax[2][1]
            corr1_ax = ax[3][0] ; corr2_ax = ax[3][1]

            sp1, fr1, bn1, im1 = spec1_ax.specgram(piezo1, Fs=1./dt, NFFT=512,
                                                   noverlap=450, xextent=xextent, mode="psd")
            sp2, fr2, bn2, im2 = spec2_ax.specgram(piezo2, Fs=1./dt, NFFT=512,
                                                   noverlap=450, xextent=xextent, mode="psd")

            spec1_ax.set_ylim(extend_window([lower, upper], 0.05))
            spec2_ax.set_ylim(extend_window([lower, upper], 0.05))

            raw1_ax.plot(timebase, piezo1)
            raw2_ax.plot(timebase, piezo2)

            n1 = len(bn1)
            n2 = len(bn2)

            sp1_sums = spectrum_sums(sp1, fr1, n1, lower, upper)
            sp2_sums = spectrum_sums(sp2, fr2, n2, lower, upper)
            sp1_sums = scipy.signal.medfilt(sp1_sums)
            sp2_sums = scipy.signal.medfilt(sp2_sums)

            rescaled_t1 = rescale_window(xextent, bn1)
            rescaled_t2 = rescale_window(xextent, bn2)

            filt1_ax.semilogy(rescaled_t1, sp1_sums)
            filt2_ax.semilogy(rescaled_t2, sp2_sums)

            corr_dt = np.mean(np.diff(rescaled_t1))
            corr_n = 1000

            corr1_signalx, corr1_signaly = corr_signal(tau, corr_dt, rescaled_t1[0], corr_n,
                                                       fit_type=piezo1_fit_type)
            corr2_signalx, corr2_signaly = corr_signal(tau, corr_dt, rescaled_t2[0], corr_n,
                                                       fit_type=piezo2_fit_type)

            corr1 = np.correlate(sp1_sums, corr1_signaly, "same")
            corr2 = np.correlate(sp2_sums, corr2_signaly, "same")

            corr1t = rescaled_t1 - 0.5 * corr_n * corr_dt
            corr2t = rescaled_t2 - 0.5 * corr_n * corr_dt
            corr1_ax.plot(corr1t, corr1, "b-")

            corr2_ax.plot(corr2t, corr2, "b-")

            inter_t01 = find_t0_from_corr(corr1t, corr1)
            inter_t02 = find_t0_from_corr(corr2t, corr2)


            # 1. Establish baseline.
            led_on = data['fastDAQ']['CAMgate'] < -0.5  # flipped definition
            first_on = np.argmax(led_on)
            first_off = np.argmin(led_on[first_on:])
            second_on = np.argmax(led_on[first_off:])
            t_first_off = timebase[first_off]
            t_second_on = timebase[second_on]

            rescaled_t_first_off1_ix = np.argmin(np.abs(rescaled_t1 - t_first_off))
            rescaled_t_first_off2_ix = np.argmin(np.abs(rescaled_t2 - t_first_off))
            rescaled_t_second_on1_ix = np.argmin(np.abs(rescaled_t1 - t_second_on))
            rescaled_t_second_on2_ix = np.argmin(np.abs(rescaled_t2 - t_second_on))

            rescaled_t_first_off1 = rescaled_t1[rescaled_t_first_off1_ix]
            rescaled_t_first_off2 = rescaled_t2[rescaled_t_first_off2_ix]
            rescaled_t_second_on1 = rescaled_t1[rescaled_t_second_on1_ix]
            rescaled_t_second_on2 = rescaled_t2[rescaled_t_second_on2_ix]

            filt1_ax.axvline(rescaled_t_first_off1, color="k")
            filt1_ax.axvline(rescaled_t_second_on1, color="k")

            filt2_ax.axvline(rescaled_t_first_off2, color="k")
            filt2_ax.axvline(rescaled_t_second_on2, color="k")

            baseline1 = np.average(sp1_sums[rescaled_t_first_off1_ix+1:rescaled_t_second_on1_ix])
            baseline2 = np.average(sp2_sums[rescaled_t_first_off2_ix+1:rescaled_t_second_on2_ix])

            baseline1rms = np.sqrt(sp1_sums[rescaled_t_first_off1_ix+1:rescaled_t_second_on1_ix].dot(sp1_sums[rescaled_t_first_off1_ix+1:rescaled_t_second_on1_ix])/
                                   sp1_sums[rescaled_t_first_off1_ix + 1:rescaled_t_second_on1_ix].size)
            baseline2rms = np.sqrt(sp2_sums[rescaled_t_first_off2_ix+1:rescaled_t_second_on2_ix].dot(sp2_sums[rescaled_t_first_off2_ix+1:rescaled_t_second_on2_ix])/
                                   sp2_sums[rescaled_t_first_off2_ix + 1:rescaled_t_second_on2_ix].size)

            filt1_ax.axhline(baseline1, color="k")
            filt2_ax.axhline(baseline2, color="k")

            inter_t01_ix = np.argmin(np.abs(rescaled_t1-inter_t01))
            inter_t02_ix = np.argmin(np.abs(rescaled_t2 - inter_t02))

            rescaled_dt1 = np.mean(np.diff(rescaled_t1))
            rescaled_dt2 = np.mean(np.diff(rescaled_t2))
            t_thresh = 100e-6 # 100 microseconds
            n_points1 = int(np.floor(t_thresh/rescaled_dt1))
            n_points2 = int(np.floor(t_thresh/rescaled_dt2))
            ix1 = 0
            while True:
                if within_baseline(sp1_sums[inter_t01_ix - n_points1 - ix1:inter_t01_ix - ix1],
                                   baseline=baseline1, rms=baseline1rms):
                    break
                ix1 += 1
                if ix1 *dt > 5*t_thresh:
                    ix1 = -1
            ix2 = 0
            while True:
                if within_baseline(sp2_sums[inter_t02_ix - n_points2 - ix2:inter_t02_ix - ix2],
                                   baseline=baseline2, rms=baseline2rms):
                    break
                ix2 += 1
                if ix2*dt > 5*t_thresh:
                    ix2=-1
            if ix1 != -1:
                new_t01 = rescaled_t1[inter_t01_ix - ix1] + rescaled_dt1 / 2
            else:
                new_t01 = -1.
            if ix2 != -1:
                new_t02 = rescaled_t2[inter_t02_ix - ix2] + rescaled_dt2 / 2
            else:
                new_t02 = -1.
            print('''
    t_first_off = {} -> {} and {}
        Indexes: {} and {}
    t_second_on = {} -> {} and {}
        Indexes: {} and {}
    Baseline = {} and {}
    Baseline RMS = {} and {}
    Lookback = {} and {}
    dt = {} and {}
    baseix = {} and {}
    Old T0 = {} and {}
    New T0 = {} and {}
'''.format(t_first_off, rescaled_t_first_off1, rescaled_t_first_off2,
           rescaled_t_first_off1_ix, rescaled_t_first_off2_ix,
           t_second_on, rescaled_t_second_on1, rescaled_t_second_on2,
           rescaled_t_second_on1_ix, rescaled_t_second_on2_ix,
           baseline1, baseline2,
           baseline1rms, baseline2rms,
           n_points1, n_points2,
           rescaled_dt1, rescaled_dt2,
           inter_t01_ix-ix1, inter_t02_ix-ix2,
           inter_t01, inter_t02,
           new_t01, new_t02))

            sub0ind1 = np.where(corr1t < 0)
            sub0ind2 = np.where(corr2t < 0)
            sub0sp1_sums = sp1_sums[sub0ind1]
            sub0sp2_sums = sp2_sums[sub0ind2]

            filt1_ax.plot(corr1_signalx + (new_t01-rescaled_t1[0]),
                          0.5*max(sub0sp1_sums) * corr1_signaly,
                          "g-", linewidth=2)
            filt2_ax.plot(corr2_signalx + (new_t02-rescaled_t2[0]),
                          0.5*max(sub0sp2_sums) * corr2_signaly, "g-", linewidth=2)
            n = len(ax.flatten())
            for subaxn in range(n):
                if subaxn % 2 == 0:
                    ax.flatten()[subaxn].axvline(x=inter_t01, color="r", linewidth=1.5)
                    ax.flatten()[subaxn].axvline(x=new_t01, color="g", linewidth=1.5)
                    pass
                else:
                    ax.flatten()[subaxn].axvline(x=inter_t02, color="r", linewidth=1.5)
                    ax.flatten()[subaxn].axvline(x=new_t02, color="g", linewidth=1.5)
                ax.flatten()[subaxn].set_xlim(xextent)
            spec1_ax.set_title("Piezo1")
            spec2_ax.set_title("Piezo2")
            spec1_ax.text(0.3, 0.3, "Spectogram", horizontalalignment="center", verticalalignment="center",
                          transform=spec1_ax.transAxes, fontsize=20)
            raw1_ax.text(0.3, 0.3, "Raw Data", horizontalalignment="center", verticalalignment="center",
                         transform=raw1_ax.transAxes, fontsize=20)
            filt1_ax.text(0.3, 0.3, "Filtered Signal", horizontalalignment="center", verticalalignment="center",
                          transform=filt1_ax.transAxes, fontsize=20)
            corr1_ax.text(0.3, 0.3, "Correlation Result", horizontalalignment="center", verticalalignment="center",
                          transform=corr1_ax.transAxes, fontsize=20)
            plt.show()
            out["bubble_t0"] = np.array([new_t01, new_t02])
            return out

        if not view_plots:  # Then don't waste time generating them
            fr1, bn1, sp1 = scipy.signal.spectrogram(piezo1, fs=1./dt, nfft=512, noverlap=450,
                                                     mode="psd", window="hanning", nperseg=512)
            fr2, bn2, sp2 = scipy.signal.spectrogram(piezo2, fs=1./dt, nfft=512, noverlap=450,
                                                     mode="psd", window="hanning", nperseg=512)

            n1 = len(bn1)
            n2 = len(bn2)

            sp1_sums = spectrum_sums(sp1, fr1, n1, lower, upper)
            sp2_sums = spectrum_sums(sp2, fr2, n2, lower, upper)
            sp1_sums = scipy.signal.medfilt(sp1_sums)
            sp2_sums = scipy.signal.medfilt(sp2_sums)

            rescaled_t1 = rescale_window(xextent, bn1)
            rescaled_t2 = rescale_window(xextent, bn2)

            corr_dt1 = np.mean(np.diff(rescaled_t1))
            corr_dt2 = np.mean(np.diff(rescaled_t2))
            corr_n1 = 1000
            corr_n2 = 1000

            corr1_signalt, corr1_signaly = corr_signal(tau, corr_dt1, rescaled_t1[0],
                                                       corr_n1, fit_type=0)
            corr2_signalt, corr2_signaly = corr_signal(tau, corr_dt2, rescaled_t2[0],
                                                       corr_n2, fit_type=0)

            corr1 = np.correlate(sp1_sums, corr1_signaly, "same")
            corr2 = np.correlate(sp2_sums, corr2_signaly, "same")

            corr1t = rescaled_t1 - 0.5 * corr_n1 * corr_dt1
            corr2t = rescaled_t2 - 0.5 * corr_n2 * corr_dt2

            inter_t01 = find_t0_from_corr(corr1t, corr1)
            inter_t02 = find_t0_from_corr(corr2t, corr2)

            # 1. Establish baseline.
            led_on = data['fastDAQ']['CAMgate'] < -0.5  # flipped definition
            first_on = np.argmax(led_on)
            first_off = np.argmin(led_on[first_on:])
            second_on = np.argmax(led_on[first_off:])
            t_first_off = timebase[first_off]
            t_second_on = timebase[second_on]

            rescaled_t_first_off1_ix = np.argmin(np.abs(rescaled_t1 - t_first_off))
            rescaled_t_first_off2_ix = np.argmin(np.abs(rescaled_t2 - t_first_off))
            rescaled_t_second_on1_ix = np.argmin(np.abs(rescaled_t1 - t_second_on))
            rescaled_t_second_on2_ix = np.argmin(np.abs(rescaled_t2 - t_second_on))

            rescaled_t_first_off1 = rescaled_t1[rescaled_t_first_off1_ix]
            rescaled_t_first_off2 = rescaled_t2[rescaled_t_first_off2_ix]
            rescaled_t_second_on1 = rescaled_t1[rescaled_t_second_on1_ix]
            rescaled_t_second_on2 = rescaled_t2[rescaled_t_second_on2_ix]

            baseline1 = np.average(sp1_sums[rescaled_t_first_off1_ix + 1:rescaled_t_second_on1_ix])
            baseline2 = np.average(sp2_sums[rescaled_t_first_off2_ix + 1:rescaled_t_second_on2_ix])

            inter_t01_ix = np.argmin(np.abs(rescaled_t1 - inter_t01))
            inter_t02_ix = np.argmin(np.abs(rescaled_t2 - inter_t02))

            rescaled_dt1 = np.mean(np.diff(rescaled_t1))
            rescaled_dt2 = np.mean(np.diff(rescaled_t2))
            t_thresh = 100e-6  # 100 microseconds
            n_points1 = int(np.floor(t_thresh / rescaled_dt1))
            n_points2 = int(np.floor(t_thresh / rescaled_dt2))

            ix1 = 0
            while True:
                if within_baseline(sp1_sums[inter_t01_ix - n_points1 - ix1:inter_t01_ix - ix1],
                                   baseline=baseline1, percent=00):
                    break
                ix1 += 1
                if ix1 *dt > 5*t_thresh:
                    ix1 = -1
            ix2 = 0
            while True:
                if within_baseline(sp2_sums[inter_t02_ix - n_points2 - ix2:inter_t02_ix - ix2],
                                   baseline=baseline2, percent=00):
                    break
                ix2 += 1
                if ix2*dt > 5*t_thresh:
                    ix2=-1
            if ix1 != -1:
                new_t01 = rescaled_t1[inter_t01_ix - ix1] + rescaled_dt1 / 2
            else:
                new_t01 = -1.
            if ix2 != -1:
                new_t02 = rescaled_t2[inter_t02_ix - ix2] + rescaled_dt2 / 2
            else:
                new_t02 = -1.

            out["bubble_t0"] = np.array([new_t01, new_t02])
            return out
    except Exception as e:
        out["bubble_t0"] = np.array([-1., -1.])
        return out



def BuildEventList(rundir, first_event=0, last_event=-1):
    # Inputs:
    #   rundir: Directory for the run
    #   first_event: Index of first event
    #   last_event: Index of last_event
    # Outputs: A sorted list of events from rundir
    # Possible Issues: If we ever move to a different naming scheme, this needs to be re-written.
    eventdirlist = os.listdir(rundir)
    eventdirlist = filter(lambda fn: (not re.search('^\d+$', fn) is None) and
                                     os.path.isdir(os.path.join(rundir, fn)),
                          eventdirlist)
    eventdirlist = filter(lambda fn: os.path.exists(os.path.join(rundir,
                                                                 *[fn, 'Event.txt'])), eventdirlist)
    eventlist = np.intp(list(eventdirlist))
    eventlist = eventlist[eventlist >= first_event]
    if last_event >= 0:
        eventlist = eventlist[eventlist <= last_event]
    return np.sort(eventlist)




def ProcessSingleRun2(rundir, dataset='SBC-2017', recondir='.', process_list=None):
    # Inputs:
    #   rundir: Location of raw data
    #   dataset: Indicator used for filtering which analyses to run
    #   recondir: Location of recon data/where we want to output our binary files
    #   process_list: List of analyses modules to run. example: ["acoustic", "event", "images"]
    # Outputs: Nothing. Saves binary files to recondir.
    if process_list is None:
        process_list = []  # This is needed since lists are mutable objects. If you have a default argument
                           # as a mutable object, then the default argument can *change* across multiple
                           # function calls since the argument is created ONCE when the function is defined.
    runname = os.path.basename(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)
    #run_recondir = os.path.join(recondir, runname)
    run_recondir = recondir

    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    if dataset == 'SBC-2017':
        loadlist = []
        for process in process_list:
            if any(process.lower().strip() == x for x in ['dytran', 'acoustic']):
                loadlist.append('fastDAQ')
        loadlist = list(set(loadlist))
    else:
        loadlist = ['~']

    acoustic_out = []
    acoustic_default = calculate_t0(None, None, None, None)


    process_list = [p.lower().strip() for p in process_list]
    print("Starting run " + rundir)
    eventlist = BuildEventList(rundir)

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + runname + '/' + str(ev))
        npev = np.array([ev], dtype=np.int32)
        thisevent = GE(rundir, ev, *loadlist)
        print('Time to load event:  '.rjust(35) +
              str(time.time() - t0) + ' seconds')


        if "acoustic" in process_list:
            # Acoustic analysis
            t1 = time.time()
            if dataset == 'SBC-2017':
                tau_peak = 0.0025884277467056165  # <-- This comes from TauResultAnalysis.py (J.G.)
                tau_average = 0.0038163479219674467  # <-- This also ^^
                lower_f = 20000
                upper_f = 40000
                piezo1_fit_type = 0
                piezo2_fit_type = 0  # See documentation of corr_t0_finder for desc of fit_types
                try:
                    acoustic_out.append(calculate_t0(data=thisevent, tau=tau_average, piezo1_fit_type=piezo1_fit_type,
                                           piezo2_fit_type=piezo2_fit_type, lower=lower_f, upper=upper_f, view_plots=False))
                except:
                    acoustic_out.append(copy.deepcopy(acoustic_default))
                acoustic_out[-1]['runid'] = runid
                acoustic_out[-1]['ev'] = npev
            et = time.time() - t1
            print('Acoustic analysis:  '.rjust(35) + str(et) + ' seconds')


        print('*** Full event analysis ***  '.rjust(35) +
              str(time.time() - t0) + ' seconds')


    if "acoustic" in process_list:
        WB(os.path.join(run_recondir,
                        'AcousticAnalysis_' + runname + '.bin'), acoustic_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    return



if __name__ == "__main__":
    import time
    from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE
    import os

    # 1. Specify directory

    # 2. Specify event_id
    # 3. Specify event number
    # 4. Compile run path
    # 5. Load the event with SBC's GetEvent function.
    # event_data = GE(ev_path, ev, "fastDAQ")
    # 6. Run! , We'll set view plots to True!
    # output = calculate_t0(data = event_data,
    #                       tau = 0.003,
    #                       lower = 20000,
    #                       upper = 40000,
    #                       view_plots = False)
    # # 7. Output looks like...
    # print("Output from run(...) is a dictionary that contains various important parameters.")
    # for k,v in output.items():
    #     print("{}: {}".format(k, v))

    runs = get_runs("/bluearc/storage/SBC-17-data/")
    runs = trim_runlist(runs, start="20171006_0", stop="20171008_0")
    for run in runs:
        if run == "20171006_3":
            evoi = 38
        elif run == "20171007_6":
            evoi = 3
        else:
            continue
        print()
        print(run)

        ev_path = os.path.join("/bluearc/storage/SBC-17-data/", run)
        data = GE(ev_path, evoi, "fastDAQ")
        calculate_t0(data, 0.0025884277, piezo1_fit_type=0, piezo2_fit_type=0, lower=20000, upper=40000, view_plots=True)
        #ProcessSingleRun2(ev_path, recondir="/pnfs/coupp/persistent/grid_output/SBC-17/T0Test/",
        #                  process_list=["acoustic"])
        # t1 = time.time()
        # event_data = GE(ev_path, str(), "fastDAQ")
        # t2 = time.time()
        # print("\tTook {}s to load event.".format(t2-t1))
        # d = calculate_t0(data = event_data,
        #                      tau = 0.0025884277,
        #                      piezo1_fit_type = 0,
        #                      piezo2_fit_type = 0,
        #                      lower = 20000,
        #                      upper = 40000,
        #                  view_plots=True)
        # t3 = time.time()
        # print("\tTook {}s to analyze entire event.".format(t3-t2))

