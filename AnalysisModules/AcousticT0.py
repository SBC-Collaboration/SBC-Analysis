from __future__ import division

import copy
import re

import numpy as np
import scipy.signal
from scipy import optimize
# import matplotlib.pyplot as plt

def my_rms(arr):
    #return np.sqrt(arr.dot(arr)/arr.size)
    return np.std(arr)

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


def calculate_t0(piezo_waveform, piezo_timebase, led_on, tau,
                 lower=20000, upper=40000, piezo_fit_type=0):
    # Inputs:
    #   piezo_waveform: A piezo waveform, generally this should have the LED pulses subtracted
    #   piezo_timebase: The times of each element in the piezo_waveform
    #   tau: The time constant we are trying to fit to the exponential decay that occurs
    #        immediately after the bubble forms
    #   lower: The lower frequency threshold for cutting off the spectrogram
    #   upper: The upper frequency threshold for cutting off the spectrogram
    #   piezo_fit_type: The type of fit to use when trying to match the filtered piezo signal. Defaults to 0.
    #   view_plots: Boolean. If true, will display some plots for analysis.
    # Outputs: A dictionary of results for the Acoustic Analysis.
    try:
        timebase = piezo_timebase
        textent = [min(timebase), max(timebase)]
        dt = np.mean(np.diff(timebase))
        fr, bn, sp = scipy.signal.spectrogram(piezo_waveform, fs=1./dt, nfft=512, noverlap=450,
                                              mode="psd", window="hanning", nperseg=512)
        n = len(bn)
        sp_sums = spectrum_sums(sp, fr, n, lower, upper)
        sp_sums = scipy.signal.medfilt(sp_sums)
        rescaled_t = rescale_window(textent, bn)
        corr_dt = np.mean(np.diff(rescaled_t))
        corr_n = 1000
        corr_t, corr_y = corr_signal(tau, corr_dt, rescaled_t[0], corr_n, fit_type=piezo_fit_type)
        corr = np.correlate(sp_sums, corr_y, "same")
        corr_t = rescaled_t - 0.5 * corr_n * corr_dt
        test_t0 = find_t0_from_corr(corr_t, corr) # This is the t0 we begin to look backwards from
        # Establish a baseline for our lookback algorithm
        #   But first we take the log of the [integrated] spectrogram signal
        log_sp_sums = np.log(sp_sums)
        first_on = np.argmax(led_on)
        first_off = np.argmin(led_on[first_on:])
        second_on = np.argmax(led_on[first_off:])
        t_first_off = timebase[first_off]
        t_second_on = timebase[second_on]
        if not np.any(led_on):
            return np.nan
        rescaled_t_first_off_index = np.argmin(np.abs(rescaled_t - t_first_off))
        rescaled_t_second_on_index = np.argmin(np.abs(rescaled_t - t_second_on))

        baseline = np.average(log_sp_sums[rescaled_t_first_off_index+1:rescaled_t_second_on_index])

        baseline_rms = my_rms(log_sp_sums[rescaled_t_first_off_index+1:rescaled_t_second_on_index])

        test_t0_index = np.argmin(np.abs(rescaled_t - test_t0))
        rescaled_dt = np.mean(np.diff(rescaled_t))
        t_thresh = 100e-6
        n_lookback = int(np.floor(t_thresh/rescaled_dt))
        pts_lookbacked_sofar = 0

        while True:
            to_test = log_sp_sums[test_t0_index-n_lookback-pts_lookbacked_sofar:test_t0_index-pts_lookbacked_sofar]

            if np.all(to_test<(baseline+5*baseline_rms)):
                break
            pts_lookbacked_sofar += 1
            if test_t0_index-n_lookback-pts_lookbacked_sofar <= 0:
                pts_lookbacked_sofar = -1
                break
        if pts_lookbacked_sofar != -1:
            t0 = rescaled_t[test_t0_index-pts_lookbacked_sofar] + rescaled_dt/2
            # plt.ioff()
            # plt.plot(rescaled_t, log_sp_sums, color="cyan")
            # plt.axhline(baseline, color="b")
            # plt.axhline(baseline+5*baseline_rms, color="k")
            # plt.axvline(t0, color="g")
            # plt.axvline(test_t0, color="r")
            # plt.axvline(t_first_off, color="k")
            # plt.axvline(t_second_on, color="k")
            # plt.show()

        else:
            t0 = np.nan
        return t0
    except Exception as e:
        return np.nan



def BandPass2(yd, f_low, f_high):
    fband = np.array([f_low, f_high])
    b, a = scipy.signal.butter(2, fband / (2.5e6 / 2.0), btype='bandpass', output='ba')
    yd_f = scipy.signal.filtfilt(b, a, yd)
    return yd_f




def CalcPiezoE(yd, td, t_wins, f_bins, t0):
    piezoE = np.zeros((t_wins.shape[0],
                       f_bins.shape[0] - 1),
                      dtype=np.float64) + np.nan

    if np.isnan(t0):
        return piezoE

    dt = td[1] - td[0]
    t_wins = t_wins + t0
    t_wins_ix = np.intp(np.round((t_wins - td[0]) / dt))
    t_wins_ix[t_wins_ix < 0] = 0
    t_wins_ix[t_wins_ix > td.shape[0]] = td.shape[0]

    for i_win in range(t_wins.shape[0]):
        this_yd = yd[t_wins_ix[i_win][0]:t_wins_ix[i_win][1]]
        if len(this_yd) < 2:
            continue
        fft_amp = np.fft.rfft(this_yd)
        fft_pow = (np.abs(fft_amp) ** 2) * dt / len(this_yd)

        df = 1 / (dt * len(this_yd))
        fd = df * (np.arange(len(fft_amp), dtype=np.float64) + 1)
        f_bins_ix = np.intp(np.round((f_bins / df) - 1))
        f_bins_ix[f_bins_ix < 0] = 0
        f_bins_ix[f_bins_ix > len(fft_amp)] = len(fft_amp)

        fft_en = fft_pow * (fd ** 2)
        for i_f in range(len(f_bins) - 1):
            piezoE[i_win, i_f] = df *\
                np.sum(fft_en[f_bins_ix[i_f]:f_bins_ix[i_f + 1]])

    return piezoE



def AcousticAnalysis(ev, tau, piezo_fit_type=0,
                     f_high=np.float64(40e3), f_low=np.float64(6e3),
                     led_amp=np.float64(-0.1), led_tau=np.float64(2e-4),
                     bs_win=np.float64([-0.15, -0.12]),
                     t0_win=np.float64([-0.12, 0]),
                     meansamp=np.intp(1e4), notbs_win=np.float64(2e-4),
                     t_wins=np.float64([[[-2e-2, -1e-2],
                                        [-1e-3, 9e-3],
                                        [-2e-4, 4e-3]],
                                        [[-2e-2, -1e-2],
                                         [-1e-3, 9e-3],
                                         [-2e-4, 4e-3]]],                                ),
                     f_bins=np.float64([[1e2, 1e3, 1e4, 1e5],
                                       [1e2, 1e3, 1e4, 1e5]]),
                     corr_lowerf=20000, corr_upperf=40000):
    # Inputs:
    #   ev: Event data (from GetEvent)
    #   tau: The expected time-constant of the exponential decay from the filtered piezo signal
    #   piezo1_fit_type: See corr_signal_types.py
    #   piezo2_fit_type: See corr_signal_Types.py
    #   f_high: High frequency used for calculating Piezo Energy
    #   f_low: Low frequency used for calculating Piezo Energy
    #   led_amp: ??
    #   led_tau: Exponential decay time constant of led-fall off (?? I think)
    #   bs_win: ??
    #   t0_win: ??
    #   meansamp: ??
    #   t_wins: ??
    #   f_bins: Frequency bins for calculating Piezo energy
    #   corr_lowerf: Lower frequency for cutting off the spectrogram of the Piezo signal
    #   corr_upperf: Higher frequency for cutting off the spectrogram of the Piezo signal
    # Outputs: A dictionary of values for the various variables we are interested in. For a list, take a look
    #          at default_output right below.
    default_output = dict(bubble_t0=np.zeros(2)-1,  # seconds
                          peak_t0=np.zeros(2)-1,  # seconds
                          piezoE=np.zeros((2,
                                           t_wins.shape[0],
                                           f_bins.shape[1] - 1),
                                          dtype=np.float64) + np.nan,
                          piezo_list=np.zeros(2)-1,
                          piezo_freq_binedges=f_bins,
                          piezo_time_windows=t_wins,
                          led_switch_off_amp=np.zeros((2, 2)),
                          led_tau=np.zeros((2, 1)),
                          led_switch_on_amp=np.zeros((2, 2)))
    out = default_output
    try:
        if not ev["fastDAQ"]["loaded"]:
            return default_output
        piezoname_list = []
        piezo_list = []
        for k in list(ev["fastDAQ"].keys()):
            m = re.search("Piezo(\d+)", k)
            if m:
                piezoname_list.append(m.group(0))
                piezo_list.append(int(m.group(1)))
        piezo_list = np.int32(piezo_list)
        ixx = np.argsort(piezo_list)
        piezo_list = piezo_list[ixx] # Sorted piezo list
        piezoname_list = [piezoname_list[ix] for ix in ixx]
        out["piezo_list"] = piezo_list
        out["bubble_t0"] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
        out["peak_t0"] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
        out["piezoE"] = np.zeros((piezo_list.shape[0], t_wins.shape[0], f_bins.shape[0]-1),
                                 dtype=np.float64) + np.nan
        out["led_switch_on_amp"] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
        out["led_switch_off_amp"] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
        out["led_tau"] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan

        fastdaq_time = ev["fastDAQ"]["time"]
        dt = ev["fastDAQ"]["caldata"]["dt"][0]

        t0_win_ix = np.intp(np.round((t0_win - fastdaq_time[0]) / dt))
        bs_win_ix = np.intp(np.round((bs_win - fastdaq_time[0]) / dt))

        k_lowpass = f_high * 2 * np.pi * dt
        k_highpass = f_high * 2 * np.pi * dt

        if "CAMgate" in ev["fastDAQ"].keys():
            led_on = ev["fastDAQ"]["CAMgate"] < -0.5
            led_switch = np.diff(np.int8(led_on))
            led_switch_on = led_switch == 1
            led_switch_off = led_switch == -1
            led_switch_on_time = fastdaq_time[:-1][led_switch_on]
            led_switch_off_time = fastdaq_time[:-1][led_switch_off]
        else:
            led_switch_on_time = np.zeros(0)
            led_switch_off_time = np.zeros(0)
    except:
        return out

    for i_piezo in range(piezo_list.shape[0]):
        try:

            raw_piezo = ev["fastDAQ"][piezoname_list[i_piezo]]
            filtered_piezo = BandPass2(raw_piezo - np.mean(raw_piezo[:meansamp]), f_low, f_high)
            led_k = 1./led_tau
            piezo_led_on = np.sum((fastdaq_time[:, None] > led_switch_on_time) *
                                  np.exp(led_k * np.minimum(led_switch_on_time -
                                                            fastdaq_time[:, None],
                                                            np.float64([0]))),
                                  axis=1)
            piezo_led_off = np.sum((fastdaq_time[:, None] > led_switch_off_time) *
                                   np.exp(led_k * np.minimum(led_switch_off_time -
                                                             fastdaq_time[:, None],
                                                             np.float64([0]))),
                                   axis=1)

            filtered_piezo_led_on = BandPass2(piezo_led_on, f_low, f_high)
            filtered_piezo_led_off = BandPass2(piezo_led_off, f_low, f_high)
            fit_func = lambda x, on_amp, off_amp: on_amp*filtered_piezo_led_on[x] + off_amp*filtered_piezo_led_off[x]
            tdata = np.arange(0, bs_win_ix[1])
            ydata = filtered_piezo[tdata]
            x0_guess = np.array([1., -1.])
            led_amp, conv = optimize.curve_fit(fit_func, tdata, ydata, x0_guess)

            yd = filtered_piezo - (led_amp[0]*filtered_piezo_led_on + led_amp[1]*filtered_piezo_led_off)
            led_only = led_amp[0] * piezo_led_on + led_amp[1] * piezo_led_off
            out["led_switch_on_amp"][i_piezo] = led_amp[0]
            out["led_switch_off_amp"][i_piezo] = led_amp[1]
            out["led_tau"][i_piezo] = led_tau
            peak_index = np.argmax(np.abs(yd[t0_win_ix[0]:t0_win_ix[1]])) + t0_win_ix[0]
            out["peak_t0"][i_piezo] = fastdaq_time[peak_index]
            t0 = calculate_t0(piezo_waveform=raw_piezo - led_only, piezo_timebase=fastdaq_time,
                              tau=tau, led_on=led_on, lower=corr_lowerf, upper=corr_upperf,
                              piezo_fit_type=piezo_fit_type)
            if abs(t0) < 10e-6:
                t0 = np.nan # Then it's essentially 0
            if abs(t0-fastdaq_time[0]) < 10e-6:
                t0 = np.nan
            out["bubble_t0"][i_piezo] = t0
            t0_index = closest_index(fastdaq_time, t0)

            out["piezoE"][i_piezo,:,:] = CalcPiezoE(raw_piezo - led_only, fastdaq_time,
                                                   t_wins, f_bins, fastdaq_time[t0_index])

        except:
            continue
    return out


if __name__ == "__main__":
    import os
    from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE
    events = {#"20171006_3": [38],
              "20171007_6": [9,10,11,12],
              }
    raw_dir = "/bluearc/storage/SBC-17-data/"
    tau = 0.0038163479219674467
    out = []
    for k,v in events.items():
        ev_path = os.path.join(raw_dir, k)
        for ev in v:
            ev_data = GE(ev_path, ev, "fastDAQ")
            out.append(AcousticAnalysis(ev_data, tau))

    for x in out:
        print(x["bubble_t0"])
