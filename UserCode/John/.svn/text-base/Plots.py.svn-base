import os
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
import scipy.signal
from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE

plt.rcParams['axes.facecolor']="white"
plt.rcParams['savefig.facecolor']="white"



def find_t0_from_corr(corrt, corry):
    # Inputs:
    #   corrt: Time-values of the correlation signal
    #   corry: Y-values of the correlation signal
    # Outputs: The time of the maximum in corry such that corrt is less than or equal to 0.
    n = np.where(corrt >= 0)
    corry[n] = 0
    return corrt[np.argmax(corry)]

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
    ycopy = deepcopy(y)
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


raw_dir = "/bluearc/storage/SBC-17-data/"
rec_dir = "/pnfs/coupp/persistent/grid_output/SBC-17/output/"

## 1 Spectrogram
runid = "20170710_0"
event = 11
t0 = RB("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output/{runid}/AcousticAnalysis_{runid}.bin".format(runid=runid))["bubble_t0"][event][1]
ev_path = os.path.join(raw_dir, runid)
cmap = "jet"
plot1 = False
if plot1:
    data = GE(ev_path, event, "fastDAQ")
    print("fastDAQ keys:", data["fastDAQ"].keys())
    p1data = data["fastDAQ"]["Piezo1"]
    p2data = data["fastDAQ"]["Piezo2"]
    t = data["fastDAQ"]["time"]
    textent = [t[0],  t[-1]]
    dt = np.mean(np.diff(t))
    nyquist_freq = 0.5/dt
    print("Nyquist frequency:", nyquist_freq)

    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True)

    sp1, fr1, bn1, im1 = ax[0].specgram(p1data, Fs=1./dt, NFFT=512,
                                           noverlap=450, xextent=textent, mode="psd",
                                        cmap=cmap)
    sp2, fr2, bn2, im2 = ax[1].specgram(p2data, Fs=1. / dt, NFFT=512,
                                           noverlap=450, xextent=textent, mode="psd",
                                        cmap=cmap)
    ax[0].set_title("Piezo1")
    ax[1].set_title("Piezo2")
    #
    # newp1d = p1data*12000 + 25*1e3
    # ax[0].plot(t, newp1d, alpha=0.4)
    #
    # newp2d = p2data*12000 + 25*1e3
    # ax[1].plot(t, newp2d, alpha=0.4)

    ax[0].set_ylabel("Frequency (MHz)")


    for subax in ax:
        subax.set_xlim(textent)
        subax.set_ylim([0, nyquist_freq])
        cur_yticks = subax.get_yticks()
        new_yticks = [ytk / 1e6 for ytk in cur_yticks]
        subax.set_yticklabels(new_yticks)
        subax.set_xlabel("Time (s)")
        subax.axvline(t0, color="k", linewidth=1.5, linestyle="--")

    #fig.colorbar(im2)
    [[x00, y00], [x01, y01]] = ax[0].get_position().get_points()
    [[x10, y10], [x11, y11]] = ax[1].get_position().get_points()
    pad = 0.01
    width = 0.02
    cbar_ax = fig.add_axes([x11 + pad, y10, width, y01 - y10])
    axcb = fig.colorbar(im2, cax=cbar_ax)
    #plt.tight_layout()
    #plt.savefig("Plots/SpecgramPlot.svg")
    #plt.savefig("Plots/SpecgramPlot.png")

    plt.show()

plot2 = False
if plot2:
    data = GE(ev_path, event, "fastDAQ")
    print("fastDAQ keys:", data["fastDAQ"].keys())
    p1data = data["fastDAQ"]["Piezo1"]
    p2data = data["fastDAQ"]["Piezo2"]
    t = data["fastDAQ"]["time"]
    textent = [t[0],  t[-1]]
    dt = np.mean(np.diff(t))
    nyquist_freq = 0.5/dt
    print("Nyquist frequency:", nyquist_freq)

    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, sharex=True)

    sp1, fr1, bn1, im1 = ax[0].specgram(p1data, Fs=1./dt, NFFT=512,
                                           noverlap=450, xextent=textent, mode="psd",
                                        cmap=cmap)
    sp2, fr2, bn2, im2 = ax[1].specgram(p2data, Fs=1. / dt, NFFT=512,
                                           noverlap=450, xextent=textent, mode="psd",
                                        cmap=cmap)
    ax[0].set_title("Piezo1")
    ax[1].set_title("Piezo2")
    ax[0].set_ylabel("Frequency (kHz)")
    # p2win = np.array([np.min(p2data), np.max(p2data)])
    # p2newwin = np.array([15*1e3, 45*1e3])
    # rs = rescale_window(p2win, p2newwin)
    # print(rs)

    #newp1d = p1data * 92000 + 25 * 1e3
    #ax[0].plot(t, newp1d, alpha=0.4)

    newp2d = p2data*7000 + 42*1e3
    ax[1].plot(t, newp2d, alpha=0.3)

    for subax in ax:
        subax.set_xlim(textent)
        subax.set_ylim([10000, 50000])
        cur_yticks = subax.get_yticks()
        new_yticks = [ytk / 1e3 for ytk in cur_yticks]
        subax.set_yticklabels(new_yticks)
        subax.set_xlabel("Time (s)")
        subax.axvline(t0, color="g", linewidth=3, linestyle="--")
        subax.axhline(15000, color="k", linewidth=1, linestyle="--")
        subax.axhline(35000, color="k", linewidth=1, linestyle="--")

    #fig.colorbar(im2)
    [[x00, y00], [x01, y01]] = ax[0].get_position().get_points()
    [[x10, y10], [x11, y11]] = ax[1].get_position().get_points()
    pad = 0.01
    width = 0.02
    cbar_ax = fig.add_axes([x11 + pad, y10, width, y01 - y10])
    axcb = fig.colorbar(im2, cax=cbar_ax)
    # ax[0].set_xlim([-0.027005171152368092, -0.0073160721318956817])
    # ax[1].set_xlim([-0.027005171152368092, -0.0073160721318956817])
    #plt.tight_layout()
    plt.savefig("Plots/SpecgramPlotZoom.svg")
    plt.savefig("Plots/SpecgramPlotZoom.png")

    plt.show()

plot3 = False
if plot3:
    lower = 20000
    upper = 40000
    data = GE(ev_path, event, "fastDAQ")
    print("fastDAQ keys:", data["fastDAQ"].keys())
    p1data = data["fastDAQ"]["Piezo1"]
    p2data = data["fastDAQ"]["Piezo2"]
    t = data["fastDAQ"]["time"]
    textent = [t[0], t[-1]]
    dt = np.mean(np.diff(t))
    nyquist_freq = 0.5 / dt
    print("Nyquist frequency:", nyquist_freq)

    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=False)

    fr1, bn1, sp1 = scipy.signal.spectrogram(p1data, fs=1. / dt, nfft=512, noverlap=450,
                                             mode="psd", window="hanning", nperseg=512)
    fr2, bn2, sp2 = scipy.signal.spectrogram(p2data, fs=1. / dt, nfft=512, noverlap=450,
                                             mode="psd", window="hanning", nperseg=512)

    n1 = len(bn1)
    n2 = len(bn2)

    sp1_sums = spectrum_sums(sp1, fr1, n1, lower, upper)
    sp2_sums = spectrum_sums(sp2, fr2, n2, lower, upper)
    sp1_sums = scipy.signal.medfilt(sp1_sums)
    sp2_sums = scipy.signal.medfilt(sp2_sums)

    rescaled_t1 = rescale_window(textent, bn1)
    rescaled_t2 = rescale_window(textent, bn2)

    #fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True)
    ax[0].plot(rescaled_t1, sp1_sums)
    ax[0].set_yticklabels([])
    ax[0].set_xlabel("Time (s)")
    ax[0].set_title("Piezo1")
    ax[0].axvline(t0, color="g", linewidth=3, linestyle="--")
    # ax[0].set_xticklabels([])
    ax[1].plot(rescaled_t2, sp2_sums)
    ax[1].set_yticklabels([])
    ax[1].set_xlabel("Time (s)")
    ax[1].set_title("Piezo2")
    ax[1].axvline(t0, color="g", linewidth=3, linestyle="--")

    ax[0].set_xlim([-0.15, 0])
    ax[1].set_xlim([-0.15, 0])
    every_nth = 2
    for subax in ax:
        for n, label in enumerate(subax.xaxis.get_ticklabels()):
            if n % every_nth != 0:
                label.set_visible(False)
    # ax[1].set_xticklabels([])
    # plt.subplots_adjust(wspace=0, hspace=0)
    #plt.savefig("Plots/SpecSums.png")
    #plt.savefig("Plots/SpecSums.svg")
    plt.show()

plot4 = False
if plot4:
    lower = 20000
    upper = 40000
    data = GE(ev_path, event, "fastDAQ")
    print("fastDAQ keys:", data["fastDAQ"].keys())
    p1data = data["fastDAQ"]["Piezo1"]
    p2data = data["fastDAQ"]["Piezo2"]
    t = data["fastDAQ"]["time"]
    textent = [t[0], t[-1]]
    dt = np.mean(np.diff(t))
    nyquist_freq = 0.5 / dt
    print("Nyquist frequency:", nyquist_freq)

    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(6,4))

    fr1, bn1, sp1 = scipy.signal.spectrogram(p1data, fs=1. / dt, nfft=512, noverlap=450,
                                             mode="psd", window="hanning", nperseg=512)
    fr2, bn2, sp2 = scipy.signal.spectrogram(p2data, fs=1. / dt, nfft=512, noverlap=450,
                                             mode="psd", window="hanning", nperseg=512)

    n1 = len(bn1)
    n2 = len(bn2)

    sp1_sums = spectrum_sums(sp1, fr1, n1, lower, upper)
    sp2_sums = spectrum_sums(sp2, fr2, n2, lower, upper)
    sp1_sums = scipy.signal.medfilt(sp1_sums)
    sp2_sums = scipy.signal.medfilt(sp2_sums)

    rescaled_t1 = rescale_window(textent, bn1)
    rescaled_t2 = rescale_window(textent, bn2)




    corr_dt = np.mean(np.diff(rescaled_t1))
    corr_n = 1000
    tau = 0.002588427
    corr1_signalx, corr1_signaly = corr_signal(tau, corr_dt, rescaled_t1[0], corr_n,
                                               fit_type=0)
    corr2_signalx, corr2_signaly = corr_signal(tau, corr_dt, rescaled_t2[0], corr_n,
                                               fit_type=0, shift=corr_n / 150)

    corr1 = np.correlate(sp1_sums, corr1_signaly, "same")
    corr2 = np.correlate(sp2_sums, corr2_signaly, "same")

    corr1t = rescaled_t1 - 0.5 * corr_n * corr_dt
    corr2t = rescaled_t2 - 0.5 * corr_n * corr_dt

    # fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True)
    ax.plot(corr1t, corr1, "b-")
    # Normalize
    mask1 = np.where(corr1t<=0)
    mask2 = np.where(corr2t<=0)
    corr1t = corr1t[mask1]
    corr2t = corr2t[mask2]
    corr1 = corr1[mask1]
    corr2 = corr2[mask2]
    corr2 = corr2*max(corr1)/max(corr2)
    ax.plot(corr2t, corr2, "g-")
    ax.set_yticklabels([])
    ax.set_xlabel("Time (s)")
    # ax[0].set_title("Piezo1")
    ax.axvline(t0, color="k", linewidth=1.5, linestyle="--")

    new_t01 = find_t0_from_corr(corr1t, corr1)
    new_t02 = find_t0_from_corr(corr2t, corr2)
    ax.axvline(new_t01, color="b", linewidth=1.5, linestyle=":")
    ax.axvline(new_t02, color="g", linewidth=1.5, linestyle="-.")
    plt.legend(["Correlation with Piezo1", "Correlation with Piezo2", "Handscan t0",
                "Piezo1 t0 from analysis", "Piezo2 t0 from analysis"],
               loc="upper left", fancybox=True, shadow=True)
    plt.xlim([min(textent), 0])
    plt.tight_layout()
    # ax[0].set_xticklabels([])
    # ax[1].plot(corr2t, corr2)
    # ax[1].set_yticklabels([])
    # ax[1].set_xlabel("Time (s)")
    # ax[1].set_title("Piezo2")
    # ax[1].axvline(t0, color="k", linewidth=1.5, linestyle="--")

    every_nth = 2
    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)
    # ax[1].set_xticklabels([])
    # plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig("Plots/CorrSample.png")
    plt.savefig("Plots/CorrSample.svg")
    plt.show()


if True:
    import gc
    runid = "20170928_0"

    events = [27, 21, 24,
              39, 43, 40,
              57, 69, 70]



    acoustic_data = RB(
        "/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output/{runid}/AcousticAnalysis_{runid}.bin".format(
            runid=runid))
    fig, ax = plt.subplots(nrows=3, ncols=3)
    ax = ax.flatten()
    for ev in range(len(events)):
        #gc.collect()
        print(ev)
        evdata = GE("/bluearc/storage/SBC-17-data/{runid}/".format(runid=runid), events[ev], "fastDAQ")
        time = evdata["fastDAQ"]["time"]
        piezo = evdata["fastDAQ"]["Piezo2"]
        ax[ev].plot(time, piezo, color="b", zorder=1)
        t0 = acoustic_data["bubble_t0"][events[ev]][1]
        if not np.isnan(t0):
            ax[ev].axvline(t0, color="m", linewidth=3, zorder=3)
    for a in ax:
        pass
        #a.set_xticklabels([])
        #a.set_yticklabels([])
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()