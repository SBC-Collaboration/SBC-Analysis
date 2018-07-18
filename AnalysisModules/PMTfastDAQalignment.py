import numpy as np
# import pdb
# execfile('readbinary.py')
# ev = ReadBlock('File location/PMTtraces.bin')
# s = ReadBlock('File location/fastDAQ_0.bin')


def PMTandFastDAQalignment(ev,
                           max_lookback=np.float64(3.0),
                           coarsening=np.intp(20),
                           force_old_style_pmtfda=False):

    default_output = dict(PMT_trigt0_sec=np.float64([-1]),
                          PMT_trigt0_frac=np.float64([-1]),
                          PMTfda_ncoinchits=np.int64(-1),
                          PMTfda_nfastDAQhits=np.int64(-1),
                          PMTfda_nPMThits=np.int64(-1),
                          PMTfda_lookback=np.float64(-1))

    if not (ev['fastDAQ']['loaded'] and
            ev['PMTtraces']['loaded'] and
            (ev['PMTtraces']['t0_frac'].shape[0] > 0) and
            (ev['PMTtraces']['t0_frac'][-1, 0] >= 0)):
        return default_output

    output_dict = default_output

    dt = ev['fastDAQ']['caldata']['dt'][0] # convert to a scalar

    # First make a -1 / 0 / 1 fastDAQ array
    fastDAQhittimes = np.float64(ev['fastDAQ']['PMTtrig'] > 2)
    if not np.any(fastDAQhittimes):
        return default_output

    if False:
        fastDAQhittimes[ev['fastDAQ']['CAMgate'] > -.5] = -1

    fastDAQhittime_risingix = np.nonzero(np.diff(fastDAQhittimes) == 1)[0]

    output_dict['PMTfda_nfastDAQhits'] =\
        np.int64(fastDAQhittime_risingix.shape[0])

    fastDAQhittimelist = dt * (.5 + np.float64(fastDAQhittime_risingix))

    pmthittimelist = (ev['PMTtraces']['t0_sec'][:, 0] -
                      ev['PMTtraces']['t0_sec'][-1, 0]) +\
        (ev['PMTtraces']['t0_frac'][:, 0] - ev['PMTtraces']['t0_frac'][-1, 0])

    pmthittimelist_late = pmthittimelist[pmthittimelist > -max_lookback]

    # Next build pmthittimes array
    pmthittimes_length = np.intp(np.floor(max_lookback / dt))
    pmthittimes = np.zeros((pmthittimes_length))
    pmthittimes[-1] = 1

    pmt_ix = -2
    while pmt_ix > -len(ev['PMTtraces']['t0_sec']):
        pmthittimes_ix = pmthittimes.size +\
            np.intp(np.round((ev['PMTtraces']['t0_sec'][pmt_ix, 0] -
                              ev['PMTtraces']['t0_sec'][-1, 0] +
                              ev['PMTtraces']['t0_frac'][pmt_ix, 0] -
                              ev['PMTtraces']['t0_frac'][-1, 0]) / dt))
        if pmthittimes_ix < 0 or pmthittimes_ix >= pmthittimes.size:
            break
        pmthittimes[pmthittimes_ix] = 1
        pmt_ix -= 1

    if not (force_old_style_pmtfda or
            ((pmthittimelist_late.shape[0] *
              fastDAQhittimelist.shape[0]) > 1e9)):

        deltat_mat = fastDAQhittimelist[:, np.newaxis] -\
            pmthittimelist_late
        deltat_list = deltat_mat.reshape(-1)

        deltat_counts = np.histogram(deltat_mat,
                                     np.arange(0, max_lookback, dt,
                                               np.float64))[0]
        lag_counts = deltat_counts[:-2] +\
            deltat_counts[1:-1] + deltat_counts[2:]

        lag_ix = np.argmax(lag_counts)

        lagwin_low = dt * np.float64(lag_ix)
        lagwin_high = 3.0 * dt + lagwin_low

        lag = np.mean(deltat_list[(deltat_list >= lagwin_low) *
                                  (deltat_list <= lagwin_high)])
        if False:
            # Could find t0 directly here, but we've only
            # looked at rising edges in fastDAQ, better to
            # do convolution
            t0 = ev['PMTtraces']['t0_sec'][-1, 0] +\
                ev['PMTtraces']['t0_frac'][-1, 0] -\
                lag.reshape(1) - ev['fastDAQ']['time'][0]
        else:
            # So we calcualte a_fine and skip the old coarsening
            # routine, proceeding to small-scale convolution
            a_fine = np.intp((max_lookback - lag) / dt) +\
                ev['fastDAQ']['time'].shape[0] -\
                np.intp(np.floor(.5 * coarsening))

    else:
        # This is the old, convolution style pmtfda
        # Now coarsen pmttrig and pmthittimes
        pmttrig_droppedsamples = np.mod(len(fastDAQhittimes),
                                        coarsening)
        pmttrig_coarse =\
            np.mean(np.reshape(fastDAQhittimes[pmttrig_droppedsamples:],
                               (-1, coarsening)), axis=1)

        # coarse_dt = dt * coarsening

        pmthittimes_droppedsamples = np.mod(len(pmthittimes),
                                            coarsening)
        pmthittimes_coarse =\
            np.mean(np.reshape(pmthittimes[pmthittimes_droppedsamples:],
                               (-1, coarsening)), axis=1)

        q_coarse = np.convolve(pmthittimes_coarse,
                               pmttrig_coarse[::-1], mode='full')
        # len(q_coarse) = len(pmthittimes_coarse) + len(pmttrig_coarse) - 1

        a = np.argmax(q_coarse)
        # a is the index of pmthittimes_coarse that's lined up with last
        # element in pmttrig_course at max overlap

        # b = a - len(pmttrig_coarse) + 1
        # b is the index of pmthittimes_coarse that's lined up with the first
        # element in pmttrig_course

        # q_fine = []

        a_fine = np.intp(a * coarsening +
                         np.floor(0.5 * coarsening) +
                         pmthittimes_droppedsamples)
        # print(a_fine)

    b_fine = a_fine - len(fastDAQhittimes) + 1

    # in all of the below, len(q_fine) = 2*coarsening + 1
    # and ix is the # of samples in pmthittimes past the end of fastDAQ

    if (a_fine + coarsening) >= len(pmthittimes):
        # slid past end of pmthittimes_coarse
        pmttrig_dropend = a_fine + coarsening - len(pmthittimes) + 1
        if (-1 - pmttrig_dropend + fastDAQhittimes.shape[0]) >= 0:
            q_fine =\
                np.convolve(pmthittimes[b_fine - coarsening:],
                            fastDAQhittimes[(-1 - pmttrig_dropend)::-1],
                            mode='valid')
        else:
            return default_output

    elif b_fine < coarsening:  # didn't slide fully onto pmthittimes_coarse
        pmttrig_dropstart = coarsening - b_fine
        if pmttrig_dropstart < (ev['fastDAQ']['PMTtrig'].shape[0] - 1):
            q_fine = np.convolve(pmthittimes[:a_fine + coarsening],
                                 fastDAQhittimes[:pmttrig_dropstart:-1],
                                 mode='valid')
        else:
            return default_output

    else:  # valid overlap in q_coarse
        q_fine = np.convolve(pmthittimes[b_fine - coarsening:
                                         a_fine + coarsening],
                             fastDAQhittimes[::-1], mode='valid')

    ix = len(pmthittimes) - (np.argmax(q_fine) + b_fine - coarsening +
                             len(ev['fastDAQ']['PMTtrig']))

    # ix = 0 if last element in pmthittimes
    # lines up with last element in ev['fastDAQ']['PMTtrig']
    t_lasthit = ix * dt + ev['fastDAQ']['time'][-1]

    t0 = ev['PMTtraces']['t0_sec'][:, 0][-1] +\
        ev['PMTtraces']['t0_frac'][:, 0][-1] -\
        t_lasthit

    output_dict['PMT_trigt0_sec'] = np.floor(t0)

    output_dict['PMT_trigt0_frac'] = t0 - np.floor(t0)

    t_rel = ev['PMTtraces']['t0_sec'][:, 0] +\
        ev['PMTtraces']['t0_frac'][:, 0] - t0
    ix_lasthitinfastDAQ = np.nonzero(t_rel <=
                                     ev['fastDAQ']['time'][-1])[0][-1]
    ix_firsthitinfastDAQ = np.nonzero(t_rel >=
                                      ev['fastDAQ']['time'][0])[0][0]

    ix_relevantfastDAQhits = np.intp(np.round((t_rel[ix_firsthitinfastDAQ:
                                                     (ix_lasthitinfastDAQ +
                                                      1)] -
                                               ev['fastDAQ']['time'][0]) /
                                              dt))

    output_dict['PMTfda_nPMThits'] = np.int64(ix_lasthitinfastDAQ -
                                              ix_firsthitinfastDAQ + 1)
    output_dict['PMTfda_ncoinchits'] = np.int64(
        np.sum((fastDAQhittimes[ix_relevantfastDAQhits] == 1) +
               (fastDAQhittimes[np.fmax(ix_relevantfastDAQhits - 1, 0)] ==
                1) +
               (fastDAQhittimes[np.fmin(ix_relevantfastDAQhits + 1,
                                        fastDAQhittimes.shape[0] - 1)] ==
                1)))
    output_dict['PMTfda_lookback'] = t_rel[-1]

    return output_dict
