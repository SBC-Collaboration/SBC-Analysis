# John Gresl 7/17/2018
# File name usages in: (change the text in these if you rename this file)
#                      corr_signal_types.py


import copy

import numpy as np
import scipy.signal
import sys

# Helper Functions


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
        out.append(np.trapz(spectrum[good_indices[0], subn], dx = np.mean(np.diff(fr))))
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

# Uncomment the below code if you want to test a new fit_type quickly to make sure the shape is what you want
# fit_type = 0
# t,y = corr_signal(1.5, 0.1, 0, 45, fit_type, 20) # <-- Uncomment to test different fit types
# plt.ioff()
# plt.plot(t,y)
# plt.show()
# 1/0 # <-- Just to stop the program here


def find_t0_from_corr(corrt, corry):
    # Inputs:
    #   corrt: Time-values of the correlation signal
    #   corry: Y-values of the correlation signal
    # Outputs: The time of the maximum in corry such that corrt is less than or equal to 0.
    n = np.where(corrt >= 0)
    corry[n] = 0
    return corrt[np.argmax(corry)]


def run(data, tau, lower, upper, piezo1_fit_type = 0, piezo2_fit_type = 4, view_plots=False):
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

    default_output = {"bubble_t0": np.zeros(2)-1.} # Default output incase the file can't be opened or something
    if not data["fastDAQ"]["loaded"]:
        return default_output
    out = default_output
    try:
        piezo1 = data["fastDAQ"]["Piezo1"]
        piezo2 = data["fastDAQ"]["Piezo2"]
        timebase = data["fastDAQ"]["time"]
        xextent = (min(timebase), max(timebase))
        dt = np.mean(np.diff(timebase))
    except Exception as e:
        print("Unable to load piezo data.")
        return out

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

            filt1_ax.plot(rescaled_t1, sp1_sums)
            filt2_ax.plot(rescaled_t2, sp2_sums)

            corr_dt = np.mean(np.diff(rescaled_t1))
            corr_n = 1000

            corr1_signalx, corr1_signaly = corr_signal(tau, corr_dt, rescaled_t1[0], corr_n, fit_type=piezo1_fit_type)
            corr2_signalx, corr2_signaly = corr_signal(tau, corr_dt, rescaled_t2[0], corr_n,
                                                       fit_type=piezo2_fit_type, shift=corr_n/150)

            corr1 = np.correlate(sp1_sums, corr1_signaly, "same")
            corr2 = np.correlate(sp2_sums, corr2_signaly, "same")

            corr1t = rescaled_t1 - 0.5 * corr_n * corr_dt
            corr2t = rescaled_t2 - 0.5 * corr_n * corr_dt
            corr1_ax.plot(corr1t, corr1, "b-")

            corr2_ax.plot(corr2t, corr2, "b-")

            new_t01 = find_t0_from_corr(corr1t, corr1)
            new_t02 = find_t0_from_corr(corr2t, corr2)

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
                    ax.flatten()[subaxn].axvline(x=new_t01, color="g", linewidth=1.5)
                    pass
                else:
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
            all_vars = locals()
            cur_size = sum([sys.getsizeof(v) for v in all_vars.values()])
            print("\t\t\tMemory Usage: {} Mb".format(cur_size*1e-6))
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
                                                       corr_n2, fit_type=4, shift=corr_n2/150)

            corr1 = np.correlate(sp1_sums, corr1_signaly, "same")
            corr2 = np.correlate(sp2_sums, corr2_signaly, "same")

            corr1t = rescaled_t1 - 0.5 * corr_n1 * corr_dt1
            corr2t = rescaled_t2 - 0.5 * corr_n2 * corr_dt2

            new_t01 = find_t0_from_corr(corr1t, corr1)
            new_t02 = find_t0_from_corr(corr2t, corr2)
            out["bubble_t0"] = np.array([new_t01, new_t02])
            return out
    except:
        out["bubble_t0"] = np.array([-1., -1.])
        return out



# EXAMPLE USAGE
if __name__ == "__main__":
    import os
    from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE
    # 1. Specify directory
    run_dir = "/bluearc/storage/SBC-17-data/"
    # 2. Specify event_id
    event_id = "20170630_0"
    # 3. Specify event number
    ev = "2"
    # 4. Compile run path
    ev_path = os.path.join(run_dir, event_id)
    # 5. Load the event with SBC's GetEvent function.
    event_data = GE(ev_path, ev, "fastDAQ")
    # 6. Run! , We'll set view plots to True!
    output = run(data = event_data,
                 tau = 0.003,
                 lower = 20000,
                 upper = 40000,
                 view_plots = False)
    # 7. Output looks like...
    print("Output from run(...) is a dictionary that contains varios important parameters.")
    for k,v in output.items():
        print("{}: {}".format(k, v))
