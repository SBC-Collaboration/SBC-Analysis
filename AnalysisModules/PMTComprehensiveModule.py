import numpy as np
# import pdb


def PMTcm_defaultout(n_traces, n_windows):
    default_output = dict(nPMThit=n_traces +
                          np.zeros([n_traces], dtype=np.int32),
                          iPMThit=1 +
                          np.arange(n_traces, dtype=np.int32),
                          PMT_t0_sec=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_t0_frac=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_t0_fastdaq=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_coinc=-1 +
                          np.zeros([n_traces], dtype=np.int16),
                          PMT_baseline=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_baserms=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_sat=-1 +
                          np.zeros([n_traces], dtype=np.int8),
                          PMT_area=np.nan +
                          np.zeros([n_traces, n_windows], dtype=np.float64),
                          PMT_area_nobs=np.nan +
                          np.zeros([n_traces, n_windows], dtype=np.float64),
                          PMT_min=np.nan +
                          np.zeros([n_traces, n_windows], dtype=np.float64),
                          PMT_max=np.nan +
                          np.zeros([n_traces, n_windows], dtype=np.float64),
                          PMT_pulse_area=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_pulse_height=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_pulse_tstart=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_pulse_tend=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_pulse_tpeak=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_pulse_t10=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_pulse_t90=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_nimtrig_ton=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_nimtrig_toff=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_ttltrig_ton=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          PMT_ttltrig_toff=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          pmt_nphe=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          pmt_npeaks=-1 +
                          np.zeros([n_traces], dtype=np.int32),
                          pmt_nclusters=-1 +
                          np.zeros([n_traces], dtype=np.int32),
                          pmt_maxpeak=np.nan +
                          np.zeros([n_traces], dtype=np.float64),
                          )
    return default_output


def PMTcm(ev,
          t_windows=np.float64([[100, 200], [0, 99]]),
          bssup_thresh=np.float64(0.002),
          base_samples=np.intp(80),
          pmtfda={'PMT_trigt0_sec': np.float64([-1]),
                  'PMT_trigt0_frac': np.float64([-1])},
          phe_width=np.float64(1.8e-9),
          phe_amp=np.float64(17e-3),
          phecount_thresh=np.float64(10),
          single_pe_thresh=np.float64(6e-3),
          single_pe_max=np.float64(2e-2),
          stop_hunting_thresh=np.float64(4e-3),
          breakout_thresh=np.float64(4e-3),
          isopeak_dt=np.float64(2e-8),
          max_tracecount=np.intp(1e5)
          ):

    n_windows = t_windows.shape[0]

    # *** First check that we have data ***
    if not (ev['PMTtraces']['loaded'] and
            (ev['PMTtraces']['t0_frac'].shape[0] > 0)
            ):
        default_output = PMTcm_defaultout(1, n_windows)
        default_output['nPMThit'][0] = -1
        default_output['iPMThit'][0] = -1
        return default_output

    # *** No pre-allocate the output dictionary
    n_traces = ev['PMTtraces']['traces'].shape[0]
    firsttrace = 0
    if n_traces > max_tracecount:
        firsttrace = n_traces - max_tracecount
        n_traces = max_tracecount

    out = PMTcm_defaultout(n_traces, n_windows)

    # *** Now fill out the timing outputs
    out['PMT_t0_sec'] = ev['PMTtraces']['t0_sec'][firsttrace:, 0]
    out['PMT_t0_frac'] = ev['PMTtraces']['t0_frac'][firsttrace:, 0]

    if ev['fastDAQ']['loaded'] and\
            ('PMT_trigt0_sec' in pmtfda) and\
            ('PMT_trigt0_frac' in pmtfda) and\
            (pmtfda['PMT_trigt0_frac'] >= 0):
        rel_t0_sec = ev['PMTtraces']['t0_sec'][firsttrace:, 0] -\
            pmtfda['PMT_trigt0_sec']
        rel_t0_frac = ev['PMTtraces']['t0_frac'][firsttrace:, 0] -\
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
            led_coinc = ev['fastDAQ']['CAMgate'] < -0.5  # LED is NIM trigger
            out['PMT_coinc'][in_fastdaq] = out['PMT_coinc'][in_fastdaq] +\
                2 * np.int16(led_coinc[fastdaq_ix[in_fastdaq]])

    # Pretty soon we'll enter the trace loop, so pre-calculate some things here
    dt = ev['PMTtraces']['dt'][0, 0]
    trig_t0 = 0
    gauss_xd = np.arange(-6, 6.001, 0.1, dtype=np.float64)
    center_ix = np.nonzero(np.abs(gauss_xd) < 1e-3)[0][0]
    phe_area = phe_width * phe_amp * np.sqrt(2 * np.pi)
    phe_width = phe_width / dt
    isopeak_dt = isopeak_dt / dt
    gauss_yd = np.exp(-np.square(gauss_xd) / (2 * phe_width * phe_width))
    gauss_ydconv = np.convolve(gauss_yd, gauss_yd)
    gauss_convnorm = np.max(gauss_ydconv)

    # Entering the trace loop
    for i_t in range(n_traces):

        # First order of business:  scale and merge the traces
        i_t_ev = i_t + firsttrace
        ls = ev['PMTtraces']['lost_samples'][i_t_ev, 0]
        if ls <= base_samples: # skip bad trace
            continue

        pmtV = ev['PMTtraces']['traces'][i_t_ev, 0, :ls] *\
            ev['PMTtraces']['v_scale'][i_t_ev, 0] +\
            ev['PMTtraces']['v_offset'][i_t_ev, 0]
        coarse_data = ev['PMTtraces']['traces'][i_t_ev, 1, :ls] *\
            ev['PMTtraces']['v_scale'][i_t_ev, 1] +\
            ev['PMTtraces']['v_offset'][i_t_ev, 1]

        satBool = (ev['PMTtraces']['traces'][i_t_ev, 0, :ls] <= -127) +\
            (ev['PMTtraces']['traces'][i_t_ev, 0, :ls] >= 126)

        np.copyto(pmtV, coarse_data, casting='same_kind', where=satBool)
        out['PMT_sat'][i_t] = np.int8(np.any(satBool))

        # Next, find the baseline and do subtraction/suppression
        this_base_samples = np.min([base_samples, np.intp(ls * .5)])

        bs = np.mean(pmtV[:this_base_samples])
        bs_rms = np.sqrt(np.var(pmtV[:this_base_samples]))
        out['PMT_baseline'][i_t] = bs
        out['PMT_baserms'][i_t] = bs_rms

        pmtV_bsub = pmtV - bs

        notbs_samples = abs(pmtV_bsub) >= bssup_thresh
        notbs_samples[1:] = notbs_samples[1:] + notbs_samples[:-1]
        notbs_samples[:-1] = notbs_samples[:-1] + notbs_samples[1:]
        pmtV_bsup = np.copy(pmtV_bsub)
        pmtV_bsup[~notbs_samples] = 0

        # Windowed calculations next -- finding areas and peaks by window
        ix_windows = np.intp(np.round((t_windows * 1e-9 - trig_t0) / dt))
        ix_windows[ix_windows < 0] = 0
        ix_windows[ix_windows > pmtV_bsub.shape[0]] = pmtV_bsub.shape[0]

        for i_w in range(n_windows):
            out['PMT_area'][i_t, i_w] = dt *\
                np.sum(pmtV_bsup[ix_windows[i_w, 0]:ix_windows[i_w, 1]])
            out['PMT_area_nobs'][i_t, i_w] = dt *\
                np.sum(pmtV_bsub[ix_windows[i_w, 0]:ix_windows[i_w, 1]])
            out['PMT_max'][i_t, i_w] = \
                np.max(pmtV_bsub[ix_windows[i_w, 0]:ix_windows[i_w, 1]])
            out['PMT_min'][i_t, i_w] = \
                np.min(pmtV_bsub[ix_windows[i_w, 0]:ix_windows[i_w, 1]])

        # Next is the single peak analysis, based on suppression threshold
        peak_ix = np.argmin(pmtV_bsub)
        peak_ix_start = np.nonzero(~notbs_samples[:peak_ix])[0]
        if peak_ix_start.shape[0] > 0:
            peak_ix_start = peak_ix_start[-1] + 1
        else:
            peak_ix_start = 0
        peak_ix_end = np.nonzero(~notbs_samples[peak_ix:])[0]
        if peak_ix_end.shape[0] > 0:
            peak_ix_end = peak_ix_end[0] + peak_ix
        else:
            peak_ix_end = pmtV_bsub.shape[0]

        out['PMT_pulse_area'][i_t] = dt *\
            np.sum(pmtV_bsub[peak_ix_start:peak_ix_end])
        out['PMT_pulse_height'][i_t] = -pmtV_bsub[peak_ix]
        out['PMT_pulse_tstart'][i_t] = trig_t0 + dt * peak_ix_start
        out['PMT_pulse_tend'][i_t] = trig_t0 + dt * (peak_ix_end - 1)
        out['PMT_pulse_tpeak'][i_t] = trig_t0 + dt * peak_ix
        pulse_frac = np.cumsum(pmtV_bsub[peak_ix_start:peak_ix_end]) /\
            np.max([out['PMT_pulse_area'][i_t], 1e-20])
        t10_ix = np.nonzero(pulse_frac >= 0.1)[0]
        t90_ix = np.nonzero(pulse_frac <= 0.9)[0]
        if t10_ix.shape[0] > 0:
            out['PMT_pulse_t10'][i_t] = trig_t0 +\
                dt * (peak_ix_start + t10_ix[0])
        else:
            out['PMT_pulse_t10'][i_t] = out['PMT_pulse_tpeak'][i_t]
        if t90_ix.shape[0] > 0:
            out['PMT_pulse_t90'][i_t] = trig_t0 +\
                dt * (peak_ix_start + t90_ix[-1])
        else:
            out['PMT_pulse_t90'][i_t] = out['PMT_pulse_tpeak'][i_t]

        # Now we do some bookkeeping for the coarse trace, in case this is LED
        nimtrig = np.nonzero(coarse_data < -0.5)[0]
        ttltrig = np.nonzero(coarse_data > 1.5)[0]
        if nimtrig.shape[0] > 0:
            out['PMT_nimtrig_ton'][i_t] = trig_t0 + dt * nimtrig[0]
            out['PMT_nimtrig_toff'][i_t] = trig_t0 + dt * nimtrig[-1]
        if ttltrig.shape[0] > 0:
            out['PMT_ttltrig_ton'][i_t] = trig_t0 + dt * ttltrig[0]
            out['PMT_ttltrig_toff'][i_t] = trig_t0 + dt * ttltrig[-1]

        # And now, phe counting *** THIS MUST BE LAST IN THE LOOP!!! ***
        if np.abs(out['PMT_area'][i_t, 0]) > (phecount_thresh * phe_area):
            out['pmt_nphe'][i_t] = out['PMT_area'][i_t, 0] / phe_area
            out['pmt_maxpeak'][i_t] = out['pmt_nphe'][i_t]
            continue

        old_xd = np.arange(1, pmtV_bsub.shape[0] + 1e-3, 1,
                           dtype=np.float64)
        new_xd = np.arange(old_xd[ix_windows[0, 0]],
                           old_xd[ix_windows[0, 1]], .1,
                           dtype=np.float64)
        newtrace = -np.interp(new_xd, old_xd, pmtV_bsub)
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
                thispeak[:thispeak_end] = p_amp *\
                    gauss_yd[-thispeak_start:]
            elif thispeak_end > thispeak.shape[0]:
                thispeak[thispeak_start:] = p_amp *\
                    gauss_yd[:(thispeak.shape[0] - thispeak_end)]
            else:
                thispeak[thispeak_start:thispeak_end] = p_amp *\
                    gauss_yd
            thistrace = thistrace - thispeak

            thispeakconv[:] = 0
            thispeakconv_start = p_ix - center_ix
            thispeakconv_end = p_ix + 3 * center_ix + 1
            if thispeakconv_start < 0:
                thispeakconv[:thispeakconv_end] = p_amp *\
                    gauss_ydconv[-thispeakconv_start:]
            elif thispeakconv_end > thispeakconv.shape[0]:
                thispeakconv[thispeakconv_start:] = p_amp *\
                    gauss_ydconv[:(thispeakconv.shape[0] - thispeakconv_end)]
            else:
                thispeakconv[thispeakconv_start:thispeakconv_end] = p_amp *\
                    gauss_ydconv
            thisconv = thisconv - thispeakconv

        if len(p_amp_list) == 0:
            out['pmt_nphe'][i_t] = 0
            out['pmt_npeaks'][i_t] = 0
            out['pmt_nclusters'][i_t] = 0
            out['pmt_maxpeak'][i_t] = 0
            continue

        p_amp_list = np.float64(p_amp_list)
        p_ix_list = np.intp(p_ix_list)
        out['pmt_maxpeak'][i_t] = p_amp_list[0] / phe_amp
        out['pmt_npeaks'][i_t] = p_amp_list.shape[0]
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
        out['pmt_nclusters'][i_t] = peak_amps.shape[0]

        peak_ix_start = p_ix_list[peak_starts]
        peak_ix_end = p_ix_list[peak_ends]

        singlephe = (peak_starts == peak_ends) * \
            (peak_amps < single_pe_max) * \
            ((np.abs(peak_ix_start - 1050) < 100) +
             (peak_amps > single_pe_thresh))

        phe_count = peak_amps / phe_amp
        phe_count[singlephe] = 1
        out['pmt_nphe'][i_t] = np.sum(phe_count)

        # pdb.set_trace()

    return out
