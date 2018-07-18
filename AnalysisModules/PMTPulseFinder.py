import numpy as np
import pdb


# Integrated PMT pulse area is in unit of V-ns across a 50 Ohm resistor
def PMTPulseFinder(ev,
                    base_samples=np.intp(80),
                    amp_gains=np.array([1, 1], dtype=np.float64),
                    threshold=np.float64(4)):

    # default output
    default_output = dict(
        ch=np.int32([-1]),
        iTrace=np.int32([-1]),
        baseline=np.float64([-1]),
        baserms=np.float64([-1]),
        pulse_area=np.float64([-1]),
        pulse_height=np.float64([-1]),
        pulse_istart=np.int32([-1]),
        pulse_iend=np.int32([-1]),
        pulse_ipeak=np.int32([-1])
    )

    if not (ev['PMTtraces']['loaded'] and
            (ev['PMTtraces']['t0_frac'].shape[0] > 0)
            ):
        return default_output

    # ADC to voltage with lost points removed
    ls = ev['PMTtraces']['lost_samples'][:, 0]  # two channels are identical
    lost_samples_min = np.min(ls[ls > 0])
    pmtV, satBool = scale(ev['PMTtraces']['v_offset'],
                          ev['PMTtraces']['v_scale'],
                          ev['PMTtraces']['traces'],
                          -127, 126, lost_samples_min, amp_gains)

    out = default_output

    if base_samples > lost_samples_min:
        base_samples = np.intp(lost_samples_min * .5)

    baseline = np.mean(pmtV[:, :base_samples], axis=1)
    baserms = np.sqrt(np.var(pmtV[:, :base_samples], axis=1))

#    out['PMT_t0_sec'] = ev['PMTtraces']['t0_sec'][:, 0]
#    out['PMT_t0_frac'] = ev['PMTtraces']['t0_frac'][:, 0]

    pmtV_bsub = pmtV - baseline[:, None]  # expand as column vectors

    # points over threshold
    overThresh = np.int32(pmtV_bsub < -threshold * baserms[:, None])

    # start and end indices of the peaks
    iPeakStart = np.nonzero(
        np.hstack(
            (overThresh[:, 0][:, None], np.diff(overThresh, n=1, axis=1))) == 1)
    iPeakEnd = np.nonzero(
        np.hstack(
            (np.diff(overThresh, n=1, axis=1), -overThresh[:, -1][:, None])) == -1)

    # Raise error if starts and ends not match
    if iPeakStart[0].shape != iPeakEnd[0].shape:
        raise ValueError("Shape mismatch!")

    # extend the peaks to the points where baseline are crossed
    # first left shift to the left boundaries
    dt = 1
    crossBase = np.zeros(iPeakStart[0].shape, dtype=bool)  # all false
    while not all(crossBase):
        # left shift the column values
        iPeakStart[1][~crossBase] -= dt
        # stop if crosses the baseline or index becomes less than 0
        iPeakStart[1][iPeakStart[1] < 0] = 0
        crossBase[np.logical_or(iPeakStart[1] == 0, pmtV_bsub[
                                iPeakStart] > 0)] = True

    # right shift
    crossBase = np.zeros(iPeakStart[0].shape, dtype=bool)  # all false
    while not all(crossBase):
        # left shift the column values
        iPeakEnd[1][~crossBase] += dt
        # stop if crosses the baseline or index becomes less than 0
        iPeakEnd[1][iPeakEnd[1] > pmtV.shape[1] - 1] = pmtV.shape[1] - 1
        crossBase[np.logical_or(iPeakEnd[1] == pmtV.shape[1] -
                                1, pmtV_bsub[iPeakEnd] > 0)] = True

    # merge overlap peaks
    nonOverlapStart = np.ones(iPeakStart[0].shape, dtype=bool)
    nonOverlapStart[1:] = ~np.logical_and(np.diff(iPeakStart[0]) == 0,
                                          np.diff(iPeakStart[1]) == 0)

    nonOverlapEnd = np.ones(iPeakEnd[0].shape, dtype=bool)
    nonOverlapEnd[1:] = ~np.logical_and(np.diff(iPeakEnd[0]) == 0,
                                        np.diff(iPeakEnd[1]) == 0)

    if any(np.logical_xor(nonOverlapStart, nonOverlapEnd)):
        raise ValueError("Numbers of peak starts and ends don't match")

    out['iTrace'] = iPeakStart[0][nonOverlapStart]
    out['baseline'] = baseline[iPeakStart[0][nonOverlapStart]]
    out['baserms'] = baserms[iPeakStart[0][nonOverlapStart]]
    out['pulse_istart'] = iPeakStart[1][nonOverlapStart]
    out['pulse_iend'] = iPeakEnd[1][nonOverlapStart]
    out['ch'] = np.zeros((np.sum(nonOverlapStart),), dtype=np.int32)
    out['ch'][satBool[iPeakStart[0][nonOverlapStart]]] = 1

    out['pulse_area'] = np.zeros((np.sum(nonOverlapStart),), dtype=np.float64)
    out['pulse_height'] = np.zeros(
        (np.sum(nonOverlapStart),), dtype=np.float64)
    out['pulse_ipeak'] = np.zeros((np.sum(nonOverlapStart),), dtype=np.int32)

    # calculate the baseline crossing points and do the integration
    # may be faster in a for loop since no need to store intermediate results
    i = 0
    for iPulse in np.nonzero(nonOverlapStart)[0]:
        # start crossing point

        x0 = np.interp(0,
                       [pmtV_bsub[(iPeakStart[0][iPulse], iPeakStart[1][iPulse])],
                        pmtV_bsub[(iPeakStart[0][iPulse], iPeakStart[1][iPulse] + 1)]],
                       [iPeakStart[1][iPulse], iPeakStart[1][iPulse] + 1])
        # end crossing point
        x1 = np.interp(0,
                       [pmtV_bsub[(iPeakEnd[0][iPulse], iPeakEnd[1][iPulse] - 1)],
                        pmtV_bsub[(iPeakEnd[0][iPulse], iPeakEnd[1][iPulse])]],
                       [iPeakEnd[1][iPulse] - 1, iPeakEnd[1][iPulse]])

        area = np.trapz(np.hstack((0,
                                   pmtV_bsub[iPeakStart[0][iPulse]][(iPeakStart[
                                       1][iPulse] + 1):iPeakEnd[1][iPulse]],
                                   0)),
                        np.hstack((x0,
                                   np.arange((iPeakStart[1][iPulse] + 1), iPeakEnd[1][iPulse]), x1)))
        iPeak = iPeakStart[1][iPulse] + np.ma.argmin(pmtV_bsub[iPeakStart[0][iPulse]][
            iPeakStart[1][iPulse]:iPeakEnd[1][iPulse]])

        height = pmtV_bsub[iPeakStart[0][iPulse]][iPeak]

        out['pulse_area'][i] = area
        out['pulse_height'][i] = height
        out['pulse_ipeak'][i] = iPeak
        i += 1

    return out


def scale(v_offset, v_scale, traces, minimum, maximum, lost_min,
          amp_gains=np.array([1, 1], dtype=np.float64)):
    fine_data = (traces[:, 0, :lost_min] * v_scale[:, 0, None] +
                 v_offset[:, 0, None]) / amp_gains[0]
    coarse_data = (traces[:, 1, :lost_min] * v_scale[:, 1, None] +
                   v_offset[:, 1, None]) / amp_gains[1]

    # only copy saturated points may not be correct if there is a delay between
    # two traces
    satBool = (traces[:, 0, :lost_min] <= minimum) +\
        (traces[:, 0, :lost_min] >= maximum)
    # np.copyto(fine_data, course_data, casting='same_kind', where=satBool)
    fine_data[np.any(satBool, axis=1)] = coarse_data[np.any(satBool, axis=1)]

    return fine_data, np.any(satBool, axis=1)


def PMTPulseFinder2(ev,
                   base_samples=np.intp(80),
                   amp_gains=np.array([1, 1], dtype=np.float64),
                   threshold=np.float64(3.5)):

    # default output
    default_output = dict(
        ch=np.int32([-1]),
        iTrace=np.int32([-1]),
        nPulse=np.int32([-1]),
        iPulse=np.int32([-1]),
        nSatADC=np.int32([-1]),
        baseline=np.float64([-1]),
        baserms=np.float64([-1]),
        pulse_area=np.float64([-1]),
        pulse_height=np.float64([-1]),
        pulse_istart=np.int32([-1]),
        pulse_iend=np.int32([-1]),
        pulse_ipeak=np.int32([-1])
    )

    if not (ev['PMTtraces']['loaded'] and
            (ev['PMTtraces']['t0_frac'].shape[0] > 0)):
        return default_output

    out = default_output

    # ADC to voltage with lost points removed
    ls = ev['PMTtraces']['lost_samples'][:, 0]  # two channels are identical
    lost_samples_min = np.min(ls[ls > 0])

    pmtV = scale2(ev['PMTtraces']['v_offset'],
                  ev['PMTtraces']['v_scale'],
                  ev['PMTtraces']['traces'],
                  lost_samples_min, amp_gains)

    if base_samples > lost_samples_min:
        base_samples = np.intp(lost_samples_min * .5)

    baseline = np.mean(pmtV[:, :, :base_samples], axis=2)
    baserms = np.sqrt(np.var(pmtV[:, :, :base_samples], axis=2))

#    out['PMT_t0_sec'] = ev['PMTtraces']['t0_sec'][:, 0]
#    out['PMT_t0_frac'] = ev['PMTtraces']['t0_frac'][:, 0]

    pmtV = pmtV - baseline[:, :, None]  # expand as column vectors

    # points and traces over threshold
    overThreshV = pmtV < -threshold * baserms[:, :, None]
    overThreshTrace = np.any(overThreshV, axis=2)
    # if no samples in a trace is over threshold, set its first sample
    # over threshold for ease of processing
    firstADC = overThreshV[:, :, 0]
    firstADC[~overThreshTrace] = True
    overThreshV[:, :, 0] = firstADC
    # convert to integer for further manipulations
    overThreshV = overThreshV.astype(int)
    # pdb.set_trace()
    # start and end indices of the peaks
    iPeakStart = np.nonzero(
        np.concatenate(
            (overThreshV[:, :, 0][:, :, None], np.diff(overThreshV, n=1, axis=2)), axis=2) == 1)
    iPeakEnd = np.nonzero(
        np.concatenate(
            (np.diff(overThreshV, n=1, axis=2), -overThreshV[:, :, -1][:, :, None]), axis=2) == -1)
    # pdb.set_trace()
    # extend the peaks to the points where the baseline are crossed
    # first left shift to the left boundaries
    dt = 1
    crossBase = np.zeros(iPeakStart[0].shape, dtype=bool)  # all false
    while not all(crossBase):
        # left shift the column values
        iPeakStart[2][~crossBase] -= dt
        # stop if crosses the baseline or index becomes less than 0
        iPeakStart[2][iPeakStart[2] < 0] = 0
        crossBase[np.logical_or(iPeakStart[2] == 0, pmtV[
                                iPeakStart] > 0)] = True

    # right shift
    crossBase = np.zeros(iPeakEnd[0].shape, dtype=bool)  # all false
    while not all(crossBase):
        # left shift the column values
        iPeakEnd[2][~crossBase] += dt
        # stop if crosses the baseline or index becomes less than 0
        iPeakEnd[2][iPeakEnd[2] > (pmtV.shape[2] - 1)] = pmtV.shape[2] - 1
        crossBase[np.logical_or(iPeakEnd[2] == pmtV.shape[2] - 1,
                                pmtV[iPeakEnd] > 0)] = True

    # merge overlap peaks
    nonOverlapStart = np.ones(iPeakStart[0].shape, dtype=bool)
    nonOverlapStart[1:] = ~np.logical_and(
        np.logical_and(np.diff(iPeakStart[0]) == 0, np.diff(iPeakStart[1]) == 0),
        np.diff(iPeakStart[2]) == 0)

    nonOverlapEnd = np.ones(iPeakEnd[0].shape, dtype=bool)
    nonOverlapEnd[1:] = ~np.logical_and(
        np.logical_and(np.diff(iPeakEnd[0]) == 0, np.diff(iPeakEnd[1]) == 0),
        np.diff(iPeakEnd[2]) == 0)

    if any(np.logical_xor(nonOverlapStart, nonOverlapEnd)):
        raise ValueError("Numbers of peak starts and ends don't match")

    # Output will include all pulses and the traces not over the threshold
    nLines = iPeakStart[0][nonOverlapStart].size
    out['ch'] = np.zeros((nLines,), dtype=np.int32) - 1
    out['iTrace'] = np.zeros((nLines,), dtype=np.int32) - 1
    out['nPulse'] = np.zeros((nLines,), dtype=np.int32) - 1
    out['iPulse'] = np.zeros((nLines,), dtype=np.int32) - 1
    out['nSatADC'] = np.zeros((nLines,), dtype=np.int32) - 1
    out['baseline'] = np.zeros((nLines,), dtype=np.float64) - 1
    out['baserms'] = np.zeros((nLines,), dtype=np.float64) - 1
    out['pulse_area'] = np.zeros((nLines,), dtype=np.float64) - 1
    out['pulse_height'] = np.zeros((nLines,), dtype=np.float64) - 1
    out['pulse_istart'] = np.zeros((nLines,), dtype=np.int32) - 1
    out['pulse_iend'] = np.zeros((nLines,), dtype=np.int32) - 1
    out['pulse_ipeak'] = np.zeros((nLines,), dtype=np.int32) - 1
    # pdb.set_trace()
    iTrace0 = 0
    ch0 = 0
    nPulse = 0
    iPulse = 0
    for iLine in range(nLines):
        iTrace = iPeakStart[0][nonOverlapStart][iLine]
        ch = iPeakStart[1][nonOverlapStart][iLine]

        out['ch'][iLine] = ch
        out['iTrace'][iLine] = iTrace
        out['baseline'][iLine] = baseline[iTrace, ch]
        out['baserms'][iLine] = baserms[iTrace, ch]

        if overThreshTrace[iTrace, ch]:  # has pulses
            if (iTrace == iTrace0) and (ch == ch0):  # same trace
                nPulse += 1
                iPulse += 1
            else:
                iTrace0 = iTrace
                ch0 = ch
                nPulse = 1
                iPulse = 1

            iStart = iPeakStart[2][nonOverlapStart][iLine]
            iEnd = iPeakEnd[2][nonOverlapEnd][iLine]

            # start crossing point
            x0 = np.interp(0, pmtV[iTrace, ch, iStart:(
                iStart + 2)], np.arange(iStart, iStart + 2))
            # end crossing point
            x1 = np.interp(0, pmtV[iTrace, ch, (iEnd - 1):(iEnd + 1)], np.arange(iEnd - 1, iEnd + 1))
            # pdb.set_trace()
            area = np.trapz(
                np.hstack((0, pmtV[iTrace, ch, (iStart + 1):iEnd], 0)),
                np.hstack((x0, np.arange(iStart + 1, iEnd), x1)))
            iPeak = iStart + np.ma.argmin(pmtV[iTrace, ch, iStart:iEnd])
            height = pmtV[iTrace, ch, iPeak]
            # pdb.set_trace()
            # Saturated ADC's
            nSatADC = np.sum(np.logical_or(ev['PMTtraces']['traces'][iTrace, ch, iStart:iEnd] < -127,
                                           (ev['PMTtraces']['traces'][iTrace, ch, iStart:iEnd] > 126)))

            out['nPulse'][(iLine - nPulse+1):(iLine+1)] = nPulse
            out['iPulse'][iLine] = iPulse
            out['nSatADC'][iLine] = nSatADC
            out['pulse_area'][iLine] = area
            out['pulse_height'][iLine] = height
            out['pulse_ipeak'][iLine] = iPeak
            out['pulse_istart'][iLine] = iStart
            out['pulse_iend'][iLine] = iEnd
        else:
            out['nPulse'][iLine] = 0
            out['iPulse'][iLine] = 0
            out['nSatADC'][iLine] = 0
            out['pulse_area'][iLine] = 0
            out['pulse_height'][iLine] = 0
            out['pulse_istart'][iLine] = 0
            out['pulse_iend'][iLine] = 0
            out['pulse_ipeak'][iLine] = 0

    return out


def scale2(v_offset, v_scale, traces, lost_min,
           amp_gains=np.array([1, 1], dtype=np.float64)):
    # dim0 - # of traces
    # dim1 - # of channels
    # dim2 - # of samples in each trace
    pmtV = (traces[:, :, :lost_min] * v_scale[:, :, None] +
            v_offset[:, :, None]) / amp_gains[None, :, None]

    return pmtV
