import numpy as np
# import ipdb
# import os
# import psutil
# from time import asctime
# @profile
def PMTpa(ev,
          t_windows=np.float64([[0, 80], [100, 300], [301, 1000]]),
          bssup_thresh=np.float64(0.002),
          base_samples=np.intp(80),
          pmtfda=None # Changed declaration to non-mutable type. JG 9/13/17
          ):
    if pmtfda is None:
        pmtfda = {'PMT_trigt0_sec': np.float64([-1]),
                  'PMT_trigt0_frac': np.float64([-1])}
    n_windows = t_windows.shape[0]
    default_output = dict(nPMThit=np.int32([-1]),
                          iPMThit=np.int32([-1]),
                          PMT_t0_sec=np.float64([-1]),
                          PMT_t0_frac=np.float64([-1]),
                          PMT_t0_fastdaq=np.float64([-1]),
                          PMT_baseline=np.float64([-1]),
                          PMT_baserms=np.float64([-1]),
                          PMT_area=np.nan +
                          np.zeros([1, n_windows], dtype=np.float64),
                          PMT_area_nobs=np.nan +
                          np.zeros([1, n_windows], dtype=np.float64),
                          PMT_min=np.nan +
                          np.zeros([1, n_windows], dtype=np.float64),
                          PMT_max=np.nan +
                          np.zeros([1, n_windows], dtype=np.float64),
                          PMT_pulse_area=np.float64([-1]),
                          PMT_pulse_height=np.float64([-1]),
                          PMT_pulse_tstart=np.float64([-1]),
                          PMT_pulse_tend=np.float64([-1]),
                          PMT_pulse_tpeak=np.float64([-1]),
                          PMT_pulse_t10=np.float64([-1]),
                          PMT_pulse_t90=np.float64([-1]),
                          PMT_coinc=np.int16([-1])
                          )
    if not (ev['PMTtraces']['loaded'] and
            (ev['PMTtraces']['t0_frac'].shape[0] > 0)
            ):
        return default_output

    ls = ev['PMTtraces']['lost_samples'][:, 0]
    lost_samples_min = np.min(ls[ls > 0])
    pmtV_bsub = scale2(ev['PMTtraces']['v_offset'],
                 ev['PMTtraces']['v_scale'],
                 ev['PMTtraces']['traces'],
                 -127, 126, lost_samples_min)

    out = default_output

    if base_samples > lost_samples_min:
        base_samples = np.intp(lost_samples_min * .5)

    out['PMT_baseline'] = np.mean(pmtV_bsub[:, :base_samples], axis=1)
    out['PMT_baserms'] = np.sqrt(np.var(pmtV_bsub[:, :base_samples], axis=1))
    out['PMT_t0_sec'] = ev['PMTtraces']['t0_sec'][:, 0]
    out['PMT_t0_frac'] = ev['PMTtraces']['t0_frac'][:, 0]

    notbs_samples = np.zeros(pmtV_bsub.shape, dtype=bool)
    for i in xrange(pmtV_bsub.shape[0]): # to save ram
        pmtV_bsub[i,:] = pmtV_bsub[i,:] - out['PMT_baseline'][i]
        notbs_samples[i,:] = abs(pmtV_bsub[i,:]) >= bssup_thresh
    notbs_samples[:, 1:] = notbs_samples[:, 1:] + notbs_samples[:, :-1]
    notbs_samples[:, :-1] = notbs_samples[:, :-1] + notbs_samples[:, 1:]
    bs_samples = ~notbs_samples
    # pmtV_bsup = np.copy(pmtV_bsub)
    # pmtV_bsup[bs_samples] = 0

    dt = ev['PMTtraces']['dt'][0, 0]
    trig_t0 = ev['PMTtraces']['t0'][0, 0]
    ix_windows = np.intp(np.round((t_windows * 1e-9 - trig_t0) / dt))
    ix_windows[ix_windows < 0] = 0
    ix_windows[ix_windows > pmtV_bsub.shape[1]] = pmtV_bsub.shape[1]

    n_traces = pmtV_bsub.shape[0]
    out['nPMThit'] = np.zeros(n_traces, dtype=np.int32) + n_traces
    out['iPMThit'] = np.arange(n_traces, dtype=np.int32) + 1
    out['PMT_area'] = np.zeros([n_traces, n_windows],
                               dtype=np.float64) + np.nan
    out['PMT_area_nobs'] = np.zeros([n_traces, n_windows],
                                    dtype=np.float64) + np.nan
    out['PMT_min'] = np.zeros([n_traces, n_windows], dtype=np.float64) + np.nan
    out['PMT_max'] = np.zeros([n_traces, n_windows], dtype=np.float64) + np.nan

    for i_win in range(n_windows):
        if (ix_windows[i_win, 1] > ix_windows[i_win, 0]) and\
                (ix_windows[i_win, 0] >= 0) and\
                (ix_windows[i_win, 1] <= pmtV_bsub.shape[1]):
            in_window = np.zeros((pmtV_bsub.shape[1],), dtype=bool)
            in_window[ix_windows[i_win, 0]:ix_windows[i_win, 1]] = True

            # out['PMT_area'][:, i_win] =\
            #     np.sum(pmtV_bsup[:, ix_windows[i_win, 0]:ix_windows[i_win, 1]],
            #            axis=1) * dt
            # out['PMT_area'][:, i_win] = \
            #     [np.sum(x[np.logical_and(i, in_window)])*dt for x, i in zip(pmtV_bsub, notbs_samples)]
            for i in xrange(0,n_traces):
                out['PMT_area'][i, i_win] = \
                        np.sum(pmtV_bsub[i, np.logical_and(notbs_samples[i,:], in_window)])*dt
            # pmtV_bsup = np.copy(pmtV_bsub[:, ix_windows[i_win, 0]:ix_windows[i_win, 1]])
            # pmtV_bsup[bs_samples[:, ix_windows[i_win, 0]:ix_windows[i_win, 1]]] = 0
            # out['PMT_area'][:, i_win] =\
            #     np.sum(pmtV_bsup,
            #            axis=1) * dt
            out['PMT_area_nobs'][:, i_win] =\
                np.sum(pmtV_bsub[:, ix_windows[i_win, 0]:ix_windows[i_win, 1]],
                       axis=1) * dt
            out['PMT_max'][:, i_win] =\
                np.max(pmtV_bsub[:, ix_windows[i_win, 0]:ix_windows[i_win, 1]],
                       axis=1)
            out['PMT_min'][:, i_win] =\
                np.min(pmtV_bsub[:, ix_windows[i_win, 0]:ix_windows[i_win, 1]],
                       axis=1)

    # for now, pulses catch last non-bs pulse in trace
    # if np.any(~bs_samples):
    #     PMT_pulse_end = np.nonzero(~bs_samples)[0][-1] + 1
    #     if np.any(bs_samples[:PMT_pulse_end]):
    #         PMT_pulse_start =\
    #             np.nonzero(bs_samples[:PMT_pulse_end])[0][-1] + 1
    #     else:
    #         PMT_pulse_start = 0
    # else:
    #     PMT_pulse_start = 0
    #     PMT_pulse_end = pmtV_bsub.shape[1]
    ixvector = np.arange(pmtV_bsub.shape[1], dtype=np.intp)

    PMT_pulse_start = np.zeros(n_traces, dtype=np.intp)
    PMT_pulse_end = np.zeros(n_traces, dtype=np.intp) + pmtV_bsub.shape[1]
    for i in xrange(0, n_traces):
        p_ends = np.nonzero((np.cumsum(~bs_samples[i, ::-1]) == 1)[::-1])
        p_starts = np.nonzero(np.cumsum(~bs_samples[i, :]) == 1)
        if p_ends[0].size:
            PMT_pulse_end[i] = p_ends[0][0] + 1
        if p_starts[0].size:
            PMT_pulse_start[i] = p_starts[0][0] + 1


    # # p_ends = np.nonzero((np.cumsum(np.cumsum(~bs_samples[:, ::-1],
    # #                                axis=1), axis=1) == 1)[:, ::-1])
    # p_ends = np.nonzero((np.cumsum(~bs_samples[:, ::-1], axis=1) == 1)[:, ::-1])
    # PMT_pulse_end[p_ends[0]] = p_ends[1] + 1
    # # p_starts = np.nonzero((np.cumsum(np.cumsum((bs_samples *
    # #                                             (ixvector <
    # #                                              (PMT_pulse_end[:, None] - 1)))
    # #                                            [:, ::-1], axis=1),
    # #                                  axis=1) == 1)[:, ::-1])
    # p_starts = np.nonzero(np.cumsum(~bs_samples, axis=1) == 1)
    # PMT_pulse_start[p_starts[0]] = p_starts[1] + 1

    # in_pulse = (ixvector >= PMT_pulse_start[:, None]) *\
    #     (ixvector < PMT_pulse_end[:, None])
    # pmtV_pulse = np.copy(pmtV_bsub)
    # pmtV_pulse[~in_pulse] = 0
    # ipdb.set_trace()
    # out['PMT_pulse_area'] = np.sum(pmtV_bsub*in_pulse, axis=1)
    # out['PMT_pulse_height'] = np.amin(pmtV_bsub, axis=1)
    out['PMT_pulse_tstart'] = dt * PMT_pulse_start
    out['PMT_pulse_tend'] = dt * PMT_pulse_end
    # out['PMT_pulse_tpeak'] = dt * np.argmin(pmtV_bsub, axis=1)
    # denom = np.copy(out['PMT_pulse_area'])
    # denom[denom == 0] = 1e-20
    # pulse_frac = np.cumsum(pmtV_bsub*in_pulse, axis=1) / denom[:, None]
    # t10s = np.nonzero(np.cumsum(pulse_frac >= 0.1, axis=1) == 1)
    # t90s = np.nonzero(np.cumsum(pulse_frac >= 0.9, axis=1) == 1)
    # pmtV_bsub = np.cumsum(pmtV_bsub*in_pulse, axis=1) / denom[:, None] # reuse pmtV_bsub to save ram
    out['PMT_pulse_height'] = np.zeros(n_traces, dtype=np.float64) + np.nan
    out['PMT_pulse_tpeak'] = np.zeros(n_traces, dtype=np.float64) + np.nan
    out['PMT_pulse_area'] = np.zeros(n_traces, dtype=np.float64)
    out['PMT_pulse_t10'] = np.zeros(n_traces, dtype=np.float64) + np.nan
    out['PMT_pulse_t90'] = np.zeros(n_traces, dtype=np.float64) + np.nan
    for i in xrange(0, n_traces):
        this_trace = pmtV_bsub[i,:]
        in_pulse = (ixvector >= PMT_pulse_start[i]) * \
                   (ixvector < PMT_pulse_end[i])
        this_trace[~in_pulse] = 0
        # out['PMT_pulse_height'][i] = np.amin(pmtV_bsub[i,in_pulse[i,:]])
        out['PMT_pulse_height'][i] = np.amin(this_trace)
        out['PMT_pulse_tpeak'][i] = np.argmin(this_trace) * dt

        # out['PMT_pulse_area'][i] = np.sum(pmtV_bsub[i,in_pulse[i,:]])
        out['PMT_pulse_area'][i] = np.sum(this_trace)
        denom = out['PMT_pulse_area'][i]
        if denom == 0:
            denom = 1e-20
        # ipdb.set_trace()
        pulse_frac = np.cumsum(this_trace)/denom
        t10s = np.nonzero(np.cumsum(pulse_frac >= 0.1) == 1)
        t90s = np.nonzero(np.cumsum(pulse_frac >= 0.9) == 1)
        if t10s[0].size:
            out['PMT_pulse_t10'][i] = dt * t10s[0][-1]
        if t90s[0].size:
            out['PMT_pulse_t90'][i] = dt * t90s[0][-1]

    # t10s = np.nonzero(np.cumsum(pmtV_bsub >= 0.1, axis=1) == 1)
    # t90s = np.nonzero(np.cumsum(pmtV_bsub >= 0.9, axis=1) == 1)

    # out['PMT_pulse_t10'] = np.zeros(n_traces, dtype=np.float64) + np.nan
    # out['PMT_pulse_t10'][t10s[0]] = dt * t10s[1]
    # out['PMT_pulse_t90'] = np.zeros(n_traces, dtype=np.float64) + np.nan
    # out['PMT_pulse_t90'][t90s[0]] = dt * t90s[1]

    # now coincidence stuff
    out['PMT_coinc'] = np.zeros(n_traces, dtype=np.int16) - 1
    out['PMT_t0_fastdaq'] = np.zeros(n_traces, dtype=np.float64)
    if ev['fastDAQ']['loaded'] and\
            ('PMT_trigt0_sec' in pmtfda) and\
            ('PMT_trigt0_frac' in pmtfda) and\
            (pmtfda['PMT_trigt0_frac'] >= 0):
        rel_t0_sec = ev['PMTtraces']['t0_sec'][:, 0] -\
            pmtfda['PMT_trigt0_sec']
        rel_t0_frac = ev['PMTtraces']['t0_frac'][:, 0] -\
            pmtfda['PMT_trigt0_frac']
        rel_t0 = rel_t0_sec + rel_t0_frac
        out['PMT_t0_fastdaq'] = rel_t0

        fastdaq_ix = np.intp(np.round((rel_t0 -
                                       ev['fastDAQ']['time'][0]) /
                                      ev['fastDAQ']['caldata']['dt']))

        in_fastdaq = (fastdaq_ix >= 0) *\
            (fastdaq_ix < ev['fastDAQ']['time'].shape[0])
        out['PMT_coinc'][in_fastdaq] = 0

        if 'VetoCoinc' in ev['fastDAQ']:
            veto_coinc = ev['fastDAQ']['VetoCoinc'] > 1.5
            veto_coinc[:-1] = veto_coinc[:-1] + veto_coinc[1:]
            veto_coinc[1:] = veto_coinc[:-1] + veto_coinc[1:]
            out['PMT_coinc'][in_fastdaq] = out['PMT_coinc'][in_fastdaq] +\
                np.int16(veto_coinc[fastdaq_ix[in_fastdaq]])

        if 'CAMgate' in ev['fastDAQ']:
            led_coinc = ev['fastDAQ']['CAMgate'] < -0.5 # LED NIM trigger negative
            out['PMT_coinc'][in_fastdaq] = out['PMT_coinc'][in_fastdaq] +\
                2 * np.int16(led_coinc[fastdaq_ix[in_fastdaq]])
    # from IPython import embed; embed()
    # pid = psutil.Process(os.getpid())
    # print(pid.memory_info().vms / 1024. ** 3)

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


# To reduce memory usage
def scale2(v_offset, v_scale, traces, minimum, maximum, lost_min):
    out_data = np.zeros((traces.shape[0], lost_min), dtype=np.float64)
    for i_trace in xrange(traces.shape[0]):
        fine_data = traces[i_trace, 0, :lost_min] * v_scale[i_trace, 0] + v_offset[i_trace, 0]
        coarse_data = traces[i_trace, 1, :lost_min] * v_scale[i_trace, 1] + v_offset[i_trace, 1]
        satBool = (traces[i_trace, 0, :lost_min] <= minimum) + \
                  (traces[i_trace, 0, :lost_min] >= maximum)
        fine_data[satBool] = coarse_data[satBool]
        out_data[i_trace, :] = fine_data
    return out_data


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
