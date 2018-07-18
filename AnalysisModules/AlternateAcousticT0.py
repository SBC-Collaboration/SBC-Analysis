from __future__ import division
import numpy as np
from scipy import signal as sig
from scipy import optimize
import re
# import ipdb


def T0finder(ev,
             f_high=np.float64(1e4), f_low=np.float64(2e3),
             led_amp=np.float64(-0.1), led_tau=np.float64(1e-4),
             bs_win=np.float64([-0.15, -0.12]),
             t0_win=np.float64([-0.12, 0]),
             meansamp=np.intp(1e4), notbs_win=np.float64(2e-4),
             t_wins=np.float64([[-2e-2, -1e-2],
                                [-1e-3, 9e-3],
                                [-2e-4, 4e-3]]),
             f_bins=np.float64([1e2, 1e3, 1e4, 1e5])):

    default_output = dict(bubble_t0=np.float64(np.nan),  # seconds
                          peak_t0=np.float64(np.nan),  # seconds
                          piezoE=np.zeros((t_wins.shape[0],
                                           f_bins.shape[0] - 1),
                                          dtype=np.float64) + np.nan,
                          piezo_list=np.int32(-1),
                          piezo_freq_binedges=f_bins,
                          piezo_time_windows=t_wins
                          )

    if not ev['fastDAQ']['loaded']:
        return default_output

    piezoname_list = []
    piezo_list = []
    for key in ev['fastDAQ'].keys():
        m = re.search('Piezo(\d+)', key)
        if m:
            piezoname_list.append(m.group(0))
            piezo_list.append(int(m.group(1)))
    piezo_list = np.int32(piezo_list)
    ixx = np.argsort(piezo_list)
    piezo_list = piezo_list[ixx]
    piezoname_list = [piezoname_list[ix] for ix in ixx]

    out = default_output
    out['piezo_list'] = piezo_list
    out['bubble_t0'] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
    out['peak_t0'] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
    out['piezoE'] = np.zeros((piezo_list.shape[0], t_wins.shape[0], f_bins.shape[0] - 1),
                             dtype=np.float64) + np.nan

    td = ev['fastDAQ']['time']
    dt = ev['fastDAQ']['caldata']['dt'][0]

    t0_win_ix = np.intp(np.round((t0_win - td[0]) / dt))
    bs_win_ix = np.intp(np.round((bs_win - td[0]) / dt))

    k_lowpass = f_high * 2 * np.pi * dt
    k_highpass = f_low * 2 * np.pi * dt

    if 'CAMgate' in ev['fastDAQ']:
        led_on = ev['fastDAQ']['CAMgate'] > -0.5
        led_switch = np.diff(np.int8(led_on))
        led_switch_on = led_switch == 1
        led_switch_off = led_switch == -1
        led_switch_on_time = td[:-1][led_switch_on]
        led_switch_off_time = td[:-1][led_switch_off]
    else:
        led_switch_on_time = np.zeros(0)
        led_switch_off_time = np.zeros(0)

    for i_piezo in range(piezo_list.shape[0]):

        ad = ev['fastDAQ'][piezoname_list[i_piezo]]

        ad_f = BandPass(ad - np.mean(ad[:meansamp]),
                        k_lowpass, k_highpass)

        led_k = 1 / led_tau
        nd = led_amp * (np.sum((td[:, None] > led_switch_on_time) *
                               np.exp(led_k * np.minimum(led_switch_on_time -
                                                         td[:, None],
                                                         np.float64([0]))),
                        axis=1) -
                        np.sum((td[:, None] > led_switch_off_time) *
                               np.exp(led_k * np.minimum(led_switch_off_time -
                                                         td[:, None],
                                                         np.float64([0]))),
                        axis=1))

        nd_f = BandPass(nd - np.mean(nd[:meansamp]),
                        k_lowpass, k_highpass)

        yd = ad_f - nd_f

        peak_ix = np.argmax(np.abs(yd[t0_win_ix[0]:t0_win_ix[1]])) + t0_win_ix[0]
        out['peak_t0'][i_piezo] = td[peak_ix]

        bs_rms = np.sqrt(np.var(yd[bs_win_ix[0]:bs_win_ix[1]]))
        not_bs = np.abs(yd) > (3 * bs_rms)
        cum_notbs = np.cumsum(np.intp(not_bs))
        notbs_winsamples = np.intp(np.round(notbs_win / dt))
        definitely_notbs = cum_notbs[notbs_winsamples:] -\
            cum_notbs[:-notbs_winsamples]
        t0_ix = np.nonzero(definitely_notbs[:(peak_ix - notbs_winsamples)] ==
                           0)[0][-1] + notbs_winsamples
        out['bubble_t0'][i_piezo] = td[t0_ix]

        out['piezoE'][i_piezo,:,:] = CalcPiezoE(ad - nd, td, t_wins, f_bins, td[t0_ix])

    return out

# T0finder2 flips the definitions of LED on and off, and calculates the LED pickup amplitudes dynamically
def T0finder2(ev,
             f_high=np.float64(40e3), f_low=np.float64(6e3),
             led_amp=np.float64(-0.1), led_tau=np.float64(2e-4),
             bs_win=np.float64([-0.15, -0.12]),
             t0_win=np.float64([-0.12, 0]),
             meansamp=np.intp(1e4), notbs_win=np.float64(2e-4),
             t_wins=np.float64([[-2e-2, -1e-2],
                                [-1e-3, 9e-3],
                                [-2e-4, 4e-3]]),
             f_bins=np.float64([1e2, 1e3, 1e4, 1e5])):

    default_output = dict(bubble_t0=np.float64(np.nan),  # seconds
                          peak_t0=np.float64(np.nan),  # seconds
                          piezoE=np.zeros((t_wins.shape[0],
                                           f_bins.shape[0] - 1),
                                          dtype=np.float64) + np.nan,
                          piezo_list=np.int32(-1),
                          piezo_freq_binedges=f_bins,
                          piezo_time_windows=t_wins
                          )

    if not ev['fastDAQ']['loaded']:
        return default_output

    piezoname_list = []
    piezo_list = []
    for key in ev['fastDAQ'].keys():
        m = re.search('Piezo(\d+)', key)
        if m:
            piezoname_list.append(m.group(0))
            piezo_list.append(int(m.group(1)))
    piezo_list = np.int32(piezo_list)
    ixx = np.argsort(piezo_list)
    piezo_list = piezo_list[ixx]
    piezoname_list = [piezoname_list[ix] for ix in ixx]

    out = default_output
    out['piezo_list'] = piezo_list
    out['bubble_t0'] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
    out['peak_t0'] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
    out['piezoE'] = np.zeros((piezo_list.shape[0], t_wins.shape[0], f_bins.shape[0] - 1),
                             dtype=np.float64) + np.nan
    out['led_switch_on_amp'] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
    out['led_switch_off_amp'] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan
    out['led_tau'] = np.zeros(piezo_list.shape, dtype=np.float64) + np.nan

    td = ev['fastDAQ']['time']
    dt = ev['fastDAQ']['caldata']['dt'][0]

    t0_win_ix = np.intp(np.round((t0_win - td[0]) / dt))
    bs_win_ix = np.intp(np.round((bs_win - td[0]) / dt))

    k_lowpass = f_high * 2 * np.pi * dt
    k_highpass = f_low * 2 * np.pi * dt

    if 'CAMgate' in ev['fastDAQ']:
        led_on = ev['fastDAQ']['CAMgate'] < -0.5 # flipped definition
        led_switch = np.diff(np.int8(led_on))
        led_switch_on = led_switch == 1
        led_switch_off = led_switch == -1
        led_switch_on_time = td[:-1][led_switch_on]
        led_switch_off_time = td[:-1][led_switch_off]
    else:
        led_switch_on_time = np.zeros(0)
        led_switch_off_time = np.zeros(0)

    for i_piezo in range(piezo_list.shape[0]):

        ad = ev['fastDAQ'][piezoname_list[i_piezo]]

        # ad_f = BandPass(ad - np.mean(ad[:meansamp]),
        #                 k_lowpass, k_highpass)
        ad_f = BandPass2(ad - np.mean(ad[:meansamp]), f_low, f_high)

        # calculate LED pickup amplitudes from baseline pickups
        pickup_len = 500
        # on_amp_min = 0.0
        # on_amp_max = 0.0
        # off_amp_min = 0.0
        # off_amp_max = 0.0
        # led_switch_on_ix = np.nonzero(led_switch_on)[0]
        # led_switch_on_ix = led_switch_on_ix[led_switch_on_ix < bs_win_ix[1]]
        # led_switch_off_ix = np.nonzero(led_switch_off)[0]

        # for jj in range(len(led_switch_on_ix)-1):
        #     on_amp_min = on_amp_min + np.amin(ad_f[led_switch_on_ix[jj] + np.arange(pickup_len)])
        #     on_amp_max = on_amp_max + np.amax(ad_f[led_switch_on_ix[jj] + np.arange(pickup_len)])
        #     off_amp_min = off_amp_min + np.amin(ad_f[led_switch_off_ix[jj] + np.arange(pickup_len)])
        #     off_amp_max = off_amp_max + np.amax(ad_f[led_switch_off_ix[jj] + np.arange(pickup_len)])
        #
        # # take the larger as amplitude
        # if np.abs(on_amp_min) > np.abs(on_amp_max):
        #     on_amp = on_amp_min / (len(led_switch_on_ix) - 1)
        # else:
        #     on_amp = on_amp_max / (len(led_switch_on_ix) - 1)
        #
        # if np.abs(off_amp_min) > np.abs(off_amp_max):
        #     off_amp = off_amp_min / (len(led_switch_on_ix) - 1)
        # else:
        #     off_amp = off_amp_max / (len(led_switch_on_ix) - 1)


        led_k = 1 / led_tau
        # nd = led_amp * (np.sum((td[:, None] > led_switch_on_time) *
        #                        np.exp(led_k * np.minimum(led_switch_on_time -
        #                                                  td[:, None],
        #                                                  np.float64([0]))),
        #                 axis=1) -
        #                 np.sum((td[:, None] > led_switch_off_time) *
        #                        np.exp(led_k * np.minimum(led_switch_off_time -
        #                                                  td[:, None],
        #                                                  np.float64([0]))),
        #                 axis=1))
        # nd = (on_amp * np.sum((td[:, None] > led_switch_on_time) *
        #                       np.exp(led_k * np.minimum(led_switch_on_time -
        #                                                 td[:, None],
        #                                                 np.float64([0]))),
        #                       axis=1) +
        #       off_amp * np.sum((td[:, None] > led_switch_off_time) *
        #                        np.exp(led_k * np.minimum(led_switch_off_time -
        #                                                  td[:, None],
        #                                                  np.float64([0]))),
        #                        axis=1))
        # nd_f = BandPass(nd - np.mean(nd[:meansamp]),
        #                 k_lowpass, k_highpass)
        # nd_f = BandPass2(nd - np.mean(nd[:meansamp]), f_low, f_high)
        # nd_f = nd_f / np.amax(nd_f) * np.abs(on_amp)

        nd_on = np.sum((td[:, None] > led_switch_on_time) *
                              np.exp(led_k * np.minimum(led_switch_on_time -
                                                        td[:, None],
                                                        np.float64([0]))),
                              axis=1)
        nd_off = np.sum((td[:, None] > led_switch_off_time) *
                               np.exp(led_k * np.minimum(led_switch_off_time -
                                                         td[:, None],
                                                         np.float64([0]))),
                               axis=1)
        nd_on_f = BandPass2(nd_on, f_low, f_high)
        nd_off_f = BandPass2(nd_off, f_low, f_high)

        func = lambda x, on_amp, off_amp: on_amp*nd_on_f[x] + off_amp*nd_off_f[x]
        xdata = np.arange(0,bs_win_ix[1])
        ydata = ad_f[xdata]
        x0 = np.array([1.0,-1.0])

        led_amp, conv = optimize.curve_fit(func, xdata, ydata, x0)
        # ipdb.set_trace()
        yd = ad_f - (led_amp[0]*nd_on_f + led_amp[1]*nd_off_f)
        nd = led_amp[0] * nd_on + led_amp[1] * nd_off

        out['led_switch_on_amp'][i_piezo] = led_amp[0]
        out['led_switch_off_amp'][i_piezo] = led_amp[1]
        out['led_tau'][i_piezo] = led_tau

        peak_ix = np.argmax(np.abs(yd[t0_win_ix[0]:t0_win_ix[1]])) + t0_win_ix[0]
        out['peak_t0'][i_piezo] = td[peak_ix]

        bs_rms = np.sqrt(np.var(yd[bs_win_ix[0]:bs_win_ix[1]]))
        not_bs = np.abs(yd) > (3.0 * bs_rms)
        cum_notbs = np.cumsum(np.intp(not_bs))
        notbs_winsamples = np.intp(np.round(notbs_win / dt))
        definitely_notbs = cum_notbs[notbs_winsamples:] -\
            cum_notbs[:-notbs_winsamples]
        t0_ix = np.nonzero(definitely_notbs[:(peak_ix - notbs_winsamples)] ==
                           0)[0][-1] + notbs_winsamples
        out['bubble_t0'][i_piezo] = td[t0_ix]

        out['piezoE'][i_piezo,:,:] = CalcPiezoE(ad - nd, td, t_wins, f_bins, td[t0_ix])

    return out

def BandPass(yd, klow, khigh):
    kl = np.exp(-klow)
    kh = np.exp(-khigh)
    yd_f = sig.lfilter([kh, -kh],
                       [1, -kh],
                       sig.lfilter([1 - kl],
                                   [1, -kl], yd))
    return yd_f


def BandPass2(yd, f_low, f_high):
    fband = np.array([f_low, f_high])
    b, a = sig.butter(2, fband / (2.5e6 / 2.0), btype = 'bandpass', output='ba')
    yd_f = sig.filtfilt(b, a, yd)
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
