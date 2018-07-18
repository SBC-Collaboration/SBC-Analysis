import numpy as np
# import pdb


def PMTphea(ev,
            t_windows=np.float64([[0, 99], [100, 200]]),
            phe_width=np.float64(1.8),
            base_samples=np.intp(80),
            single_pe_thresh=np.float64(6e-3),
            single_pe_max=np.float64(1e-3),
            stop_hunting_thresh=np.float64(2e-3),
            breakout_thresh=np.float64(2e-3),
            isopeak_dt=np.float64(20),
            phe_amp=np.float64(18e-3)
            ):
    n_windows = t_windows.shape[0]
    default_output = dict(nPMThit=np.int32([-1]),
                          iPMThit=np.int32([-1]),
                          pmt_nphe=np.nan +
                          np.zeros([1, n_windows], dtype=np.float64),
                          pmt_npeaks=-1 +
                          np.zeros([1, n_windows], dtype=np.int32),
                          pmt_nclusters=-1 +
                          np.zeros([1, n_windows], dtype=np.int32),
                          pmt_maxpeak=np.nan + np.float64([0])
                          )
    if not (ev['PMTtraces']['loaded'] and
            (ev['PMTtraces']['t0_frac'].shape[0] > 0)
            ):
        return default_output

    ls = ev['PMTtraces']['lost_samples'][:, 0]
    lost_samples_min = np.min(ls[ls > 0])
    pmtV = scale(ev['PMTtraces']['v_offset'],
                 ev['PMTtraces']['v_scale'],
                 ev['PMTtraces']['traces'],
                 -127, 126, lost_samples_min)

    out = default_output

    if base_samples > lost_samples_min:
        base_samples = np.intp(lost_samples_min * .5)

    bs = np.mean(pmtV[:, :base_samples], axis=1)

    pmtV_bsub = pmtV - bs[:, None]

    dt = ev['PMTtraces']['dt'][0, 0]
    trig_t0 = ev['PMTtraces']['t0'][0, 0]
    ix_windows = np.intp(np.round((t_windows * 1e-9 - trig_t0) / dt))
    ix_windows[ix_windows < 0] = 0
    ix_windows[ix_windows > pmtV_bsub.shape[1]] = pmtV_bsub.shape[1]
    firstsampleused = np.min(ix_windows[:, 0])
    lastsampleused = np.max(ix_windows[:, 1])

    n_traces = pmtV.shape[0]
    out['nPMThit'] = np.zeros(n_traces, dtype=np.int32) + n_traces
    out['iPMThit'] = np.arange(n_traces, dtype=np.int32) + 1
    out['pmt_nphe'] = np.zeros([n_traces, n_windows],
                               dtype=np.float64) + np.nan
    out['pmt_npeaks'] = np.zeros([n_traces, n_windows],
                                 dtype=np.int32) - 1
    out['pmt_nclusters'] = np.zeros([n_traces, n_windows],
                                    dtype=np.int32) - 1
    out['pmt_maxpeak'] = np.zeros(n_traces,
                                  dtype=np.float64) + np.nan
    gauss_xd = np.arange(-6, 6.001, 0.1, dtype=np.float64)
    center_ix = np.nonzero(np.abs(gauss_xd) < 1e-3)[0][0]
    gauss_yd = np.exp(-np.square(gauss_xd) / (2 * phe_width * phe_width))
    gauss_ydconv = np.convolve(gauss_yd, gauss_yd)
    gauss_convnorm = np.max(gauss_ydconv)

    old_xd = np.arange(1, pmtV.shape[1] + 1e-3, 1, dtype=np.float64)
    new_xd = np.arange(old_xd[firstsampleused], old_xd[lastsampleused],
                       .1, dtype=np.float64)

    for i_t in range(pmtV.shape[0]):
        newtrace = -np.interp(new_xd, old_xd, pmtV_bsub[i_t])
        thistrace = np.copy(newtrace)
        thisconv = np.convolve(thistrace, gauss_yd)
        thispeak = np.zeros(thistrace.shape, dtype=np.float64)
        thispeakconv = np.zeros(thisconv.shape, dtype=np.float64)
        p_amp_list = []
        p_ix_list = []
        while np.max(thistrace) > stop_hunting_thresh:
            p_ix = np.argmax(thisconv)
            p_amp = thisconv[p_ix]
            p_amp = p_amp / gauss_convnorm
            p_ix = p_ix - center_ix

            if p_amp < breakout_thresh or \
                    p_ix < center_ix or \
                    p_ix > (thistrace.shape[0] - center_ix - 1):
                break

            alt_traceslice = thistrace[(p_ix - center_ix):
                                       (p_ix - center_ix + gauss_yd.shape[0])]
            alternate_p_ix = np.argmax(alt_traceslice)
            alternate_max_p_amp = alt_traceslice[alternate_p_ix] * \
                np.exp(0.125 / phe_width)
            alternate_p_ix = alternate_p_ix + p_ix - center_ix

            if p_amp > alternate_max_p_amp:
                p_amp = alternate_max_p_amp
                p_ix = alternate_p_ix

            p_amp_list.append(p_amp)
            p_ix_list.append(p_ix)

            thispeak[:] = 0
            thispeak_start = p_ix - center_ix
            thispeak_end = p_ix + center_ix + 1
            if thispeak_start < 0:
                thispeak[:thispeak_end] = p_amp * gauss_yd[-thispeak_start:]
            elif thispeak_end > thispeak.shape[0]:
                thispeak[thispeak_start:] = p_amp * gauss_yd[:(thispeak.shape[0] - thispeak_end)]
            else:
                thispeak[thispeak_start:thispeak_end] = p_amp * gauss_yd
            thistrace = thistrace - thispeak

            thispeakconv[:] = 0
            thispeakconv_start = p_ix - center_ix
            thispeakconv_end = p_ix + 3 * center_ix + 1
            if thispeakconv_start < 0:
                thispeakconv[:thispeakconv_end] = p_amp * gauss_ydconv[-thispeakconv_start:]
            elif thispeakconv_end > thispeakconv.shape[0]:
                thispeakconv[thispeakconv_start:] = p_amp * gauss_ydconv[:(thispeakconv.shape[0] - thispeakconv_end)]
            else:
                thispeakconv[thispeakconv_start:thispeakconv_end] = p_amp * gauss_ydconv
            thisconv = thisconv - thispeakconv

        if len(p_amp_list) == 0:
            out['pmt_nphe'][i_t, :] = 0
            out['pmt_npeaks'][i_t, :] = 0
            out['pmt_nclusters'][i_t, :] = 0
            out['pmt_maxpeak'][i_t] = 0
            continue

        p_amp_list = np.float64(p_amp_list)
        p_ix_list = np.intp(p_ix_list)
        out['pmt_maxpeak'][i_t] = p_amp_list[0]
        timeorder = np.argsort(p_ix_list)
        p_ix_list = p_ix_list[timeorder]
        p_amp_list = p_amp_list[timeorder]

        peak_tdiff = np.concatenate((np.float64([np.inf]),
                                     np.diff(p_ix_list),
                                     np.float64([np.inf])))
        peak_starts = np.nonzero(peak_tdiff[:-1] > 10 * isopeak_dt)[0]
        peak_ends = np.nonzero(peak_tdiff[1:] > 10 * isopeak_dt)[0]

        cumsum_amps = np.zeros(peak_tdiff.shape[0], dtype=np.float64)
        cumsum_amps[1:] = np.cumsum(p_amp_list)
        cumsum_peaks = np.zeros(peak_starts.shape[0] + 1, dtype=np.float64)
        cumsum_peaks[:-1] = cumsum_amps[peak_starts]
        cumsum_peaks[-1] = cumsum_amps[-1]
        peak_amps = np.diff(cumsum_peaks)

        peak_ix_start = p_ix_list[peak_starts]
        peak_ix_end = p_ix_list[peak_ends]

        singlephe = (peak_starts == peak_ends) * \
            (peak_amps < single_pe_max) * \
            ((np.abs(peak_ix_start - 1050) < 100) +
             (peak_amps > single_pe_thresh))

        phe_count = peak_amps / phe_amp
        phe_count[singlephe] = 1

        for i_win in range(n_windows):
            thiswin_cut = ~((peak_ix_end < (10 * (ix_windows[i_win][0] - firstsampleused))) +
                            (peak_ix_start > (10 * (ix_windows[i_win][1] - firstsampleused))))
            out['pmt_nphe'][i_t, i_win] = np.sum(phe_count[thiswin_cut])
            out['pmt_nclusters'][i_t, i_win] = np.sum(thiswin_cut)
            out['pmt_npeaks'][i_t, i_win] = \
                np.sum((p_ix_list > (10 * ix_windows[i_win][0])) *
                       (p_ix_list < (10 * ix_windows[i_win][1])))
            # pdb.set_trace()
    return out


# By Allison Grimsted
# June 2016

'''
Roadmap:
    -main
    -Functions called by main:
        -scale
        -make_dictinary
    -Functions called by make_dictionary
        -integrate
        -half_integral_time
        -integral_end_time
        -integral_start_time
        -total_time
        -max_time
        -mean
        -std_dev
    -Called by secondary functions:
        -low_lim_minimum
        -up_lim_minimum
        -minimum index
        -time
'''

# import matplotlib.pyplot as plt
# from SBCcode.DataHandling.ReadBinary import ReadBlock as rb


# loads binary file
# print "Loading PMT Traces"
# dic=rb('PMTtraces20160622_0_4.bin')
# print "done loading PMT traces"


def main(event):
    default_output = {}
    npnan = np.nan + np.zeros(1, dtype=np.float64)

    default_output.update({"PMTpulse_baseline": npnan})
    default_output.update({"PMTpulse_baserms": npnan})
    default_output.update({"PMTpulse_area_specific": npnan})
    default_output.update({"PMTpulse_area_total": npnan})
    default_output.update({"PMTpulse_height": npnan})
    default_output.update({"PMTpulse_half_time_specif": npnan})
    default_output.update({"PMTpulse_half_time_tot": npnan})
    default_output.update({"PMTpulse_integral_end_time_speci": npnan})
    default_output.update({"PMTpulse_tot_time": npnan})
    default_output.update({"PMTpulse_max_voltage_time": npnan})

    if not event['PMTtraces']['loaded']:
        return default_output

    else:
        dic = event['PMTtraces']
        lost_samples = dic['lost_samples'][:, 0]
        base_vals = 80
        minimum = -127
        maximum = 126

        # finds the global minimum for the number of good samples.
        lost_min = np.min(lost_samples[lost_samples > 0])

        scaled = scale(dic['v_offset'], dic['v_scale'],
                       dic['traces'], minimum, maximum, lost_min)
        d = make_dictionary(dic['dt'][0][0], base_vals,
                            lost_min, maximum, scaled)

        # for graphing:
        '''
        dt=dic['dt'][0][0]
        x=time(dt, lost_min)
        plt.plot(x, scaled[1])
        plt.show()
        '''

        return d


def scale(v_offset, v_scale, traces, minimum, maximum, lost_min):
    fine_data = traces[:, 0, :lost_min] * v_scale[:, 0, None] +\
        v_offset[:, 0, None]
    course_data = traces[:, 1, :lost_min] * v_scale[:, 1, None] +\
        v_offset[:, 1, None]

    satBool = (traces[:, 0, :lost_min] <= minimum) +\
        (traces[:, 0, :lost_min] >= maximum)
    np.copyto(fine_data, course_data, casting='same_kind', where=satBool)

    return fine_data


def make_dictionary(dt, base_vals, lost_min, max_val, arr):
    # creates a dictionary for values to be written to a file.

    d = {}

    # declares arrays to be associated with a key.
    average = np.zeros(len(arr))
    stand = np.zeros(len(arr))
    area_specific = np.zeros(len(arr))
    area_total = np.zeros(len(arr))
    height = np.zeros(len(arr))
    half_time_specif = np.zeros(len(arr))
    half_time_tot = np.zeros(len(arr))
    integral_end_time_speci = np.zeros(len(arr))
    tot_time = np.zeros(len(arr))
    max_t = np.zeros(len(arr))

    average = np.mean(arr[:, :base_vals], axis=1)
    stand = np.std(arr[:, :base_vals], axis=1)

    for i in range(len(arr)):
        avg = average[i]
        std = stand[i]

        # fills arrays
        # average[i] = mean(base_vals, arr[i])
        # stand[i] = std_dev(base_vals, arr[i])

        low_ix = low_lim_minimum((avg - 3 * std), arr[i])
        if low_ix > 5:
            low_ix -= 5
        else:
            low_ix = 0

        high_ix = up_lim_minimum(avg, arr[i])

        area_specific[i] = integrate(dt, avg,
                                     low_ix,
                                     high_ix, arr[i])
        area_total[i] = integrate(dt, avg, 0, len(arr[i]) - 1, arr[i])

        height[i] = (avg - np.min(arr[i]))
        half_time_specif[i] =\
            half_integral_time(dt, avg,
                               low_ix,
                               high_ix, arr[i])
        half_time_tot[i] = half_integral_time(dt, avg, 1,
                                              (len(arr[i]) - 2), arr[i])
        integral_end_time_speci[i] =\
            integral_end_time(dt, avg,
                              high_ix,
                              arr[i])
        tot_time[i] = total_time(dt, lost_min)
        max_t[i] = max_time(dt, max_val, arr[i])

    # fills dictionary
    d.update({"PMTpulse_baseline": average})
    d.update({"PMTpulse_baserms": stand})
    d.update({"PMTpulse_area_specific": area_specific})
    d.update({"PMTpulse_area_total": area_total})
    d.update({"PMTpulse_height": height})
    d.update({"PMTpulse_half_time_specif": half_time_specif})
    d.update({"PMTpulse_half_time_tot": half_time_tot})
    d.update({"PMTpulse_integral_end_time_speci": integral_end_time_speci})
    d.update({"PMTpulse_tot_time": tot_time})
    d.update({"PMTpulse_max_voltage_time": max_t})

    return d


def integrate(dt, base_line, low_lim, up_lim, y):
    ''' this function assumes that the derivative of the
        line between two values of dt is constant.

        calculates the the area between the values at the
        baseline closest to the largest minimum.
        '''

    # determines the area of the right triangle between the baseline and
    # the index of the closest value to the baseline near the lower limit.
    base = base_line - y[low_lim]
    height = 0
    if (y[low_lim] - y[low_lim - 1]) != 0:
        height = (float(y[low_lim] - base_line)) /\
            (y[low_lim] - y[low_lim - 1]) * dt
    area_low = 0.5 * base * height

    # finds the area of the right triangle between the baseline and the index
    # of the closest value to the base line near the upper limit.
    area_up = 0
    if up_lim != len(y) - 1:
        h = 0
        if (y[up_lim + 1] - y[up_lim]) != 0:
            h = (base_line - y[up_lim]) * dt /\
                (y[up_lim + 1] - y[up_lim])
        b = base_line - y[up_lim]
        area_up = 0.5 * h * b

    # finds the area between the two limits found.
    total = 0
    for i in range(low_lim, up_lim):
        total += abs(0.5 * dt * ((base_line - y[i]) +
                     (base_line - y[i + 1])))
    summa = total + area_low + area_up

    return summa


def half_integral_time(dt, base_line, low_lim, up_lim, y):
    # finds the time when 50% of the light (the integral) has come through
    val = 0.5 * integrate(dt, base_line, low_lim, up_lim, y)
    t_sum = 0
    a_left = val
    i = low_lim
    t_sum = dt * low_lim

    # determines the area of the right triangle between the baseline and the
    # index of the closest value to the baseline near the lower limit.
    base = base_line - y[low_lim]
    height = 0
    if (y[low_lim] - y[low_lim - 1]) != 0:
        height = (float(y[low_lim] - base_line)) /\
            (y[low_lim] - y[low_lim - 1]) * dt
    area_low = 0.5 * base * height

    a_left = a_left - area_low

    while i < up_lim:
        area = 0.5 * dt * ((base_line - y[i]) + (base_line - y[i + 1]))
        a_left = a_left - 0.5 * dt * ((base_line - y[i]) +
                                      (base_line - y[i + 1]))
        if a_left < 0:
            a_left = a_left + 0.5 * dt * ((base_line - y[i]) +
                                          (base_line - y[i + 1]))
            break
        else:
            i += 1

    if (i - 1 - low_lim) > (-1):
        t_sum += (i - 1 - low_lim) * dt

    # calculates time for smaller piece of area.
    area_needed = a_left
    radicand = 0
    h1 = 0
    h2 = 0
    h = h2
    if (i + 1) < len(y):
        radicand = float((2 * dt * y[i] - 2 * base_line * dt) ** 2 -
                         8 * area_needed * dt * (y[i + 1] - y[i]))
        square_root = float(abs(radicand) ** 0.5)
        if (2 * y[i + 1] - 2 * y[i]) != 0:
            h1 = float(-1 * square_root + 2 * base_line * dt -
                       2 * dt * y[i]) /\
                (2 * y[i + 1] - 2 * y[i])
            h2 = float(square_root + 2 * base_line * dt - 2 * dt * y[i]) /\
                (2 * y[i + 1] - 2 * y[i])
        if h2 > dt or h2 < 0:
            h = h1
    t_sum += h
    return t_sum


def integral_end_time(dt, base_line, up_lim, y):
    # finds the time when all of the light has come through.
    time = dt * up_lim
    if up_lim != len(y) - 1:
        time += (base_line - y[up_lim]) * dt /\
            (y[up_lim + 1] - y[up_lim])
    return time


def integral_start_time(dt, low_lim):
    # finds the time when the light starts coming in.
    time = low_lim * dt
    return time


def total_time(dt, lost_min):
    # finds the time of the whole pulse
    time = dt * lost_min
    return time


def max_time(dt, max_val, y):
    # returns the time when the voltage spikes.
    # if there is not a value over the max, the function returns -1.
    time = -1
    ix = np.argmax(y)
    if y[ix] > max_val:
        time = dt * ix
    return time


'''
def mean(base_vals, y):
    # returns the mean value of the base line values

    values = y[:base_vals]

    mean = np.mean(values)
    return mean


def std_dev(base_vals, y):
    # returns the standard deviation of the base line values.

    values = y[:base_vals]

    sigma = np.std(values)
    return sigma
'''


def low_lim_minimum(base_line, y):
    '''finds the index of the lower limit by finding
    the first value to the left of the minimum
    that returns to the baseline'''

    index = np.argmin(y)
    low_lim = 0
    j = index
    caught = False
    while (not caught) and j > 1:
        if base_line == y[j - 1]:
            caught = True
            low_lim = j
        elif base_line < y[j - 1] and base_line > y[j]:
            low_lim = j
            caught = True
        else:
            j -= 1
    return low_lim


def up_lim_minimum(base_line, y):
    index = np.argmin(y)
    '''finds the index of the upper limit by finding the first value
    to the right of the minimum
    that returns to the baseline'''
    up_lim = 0
    k = index
    found = False
    if y[len(y) - 1] < base_line:
        found = True
        up_lim = len(y) - 1
    while (not found) and k < (len(y) - 1):
        if base_line == y[k + 1]:
            found = True
            up_lim = k
        elif base_line < y[k + 1] and base_line > y[k]:
            up_lim = k
            found = True
        else:
            k += 1
    return up_lim


'''
def minimum_index(y):
    # finds the index of the minimum
    if True:
        return np.argmin(y)

    minimum = np.min(y)
    i = 0
    index = 0
    found = False
    while (not found) and i < len(y):
        if y[i] == minimum:
            found = True
            index = i
        i += 1
    return index
'''


def time(dt, maximum):
    # time_elapsed is a list that contains time values based on dt
    '''time_elapsed = []
    for i in range(maximum):
        time_elapsed.append(i * dt)'''
    time_elapsed = np.arange(maximum) * dt

    return time_elapsed
