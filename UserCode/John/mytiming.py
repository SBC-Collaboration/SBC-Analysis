import numpy as np
# import ipdb

def TimingAnalysis(ev, pmtpa, aa):
    n_piezos = 2
    n_windows = 2
    default_output = dict(CAMperiod=np.float64(-1),  # seconds
                          CAMduty=np.float64(-1),  # seconds
                          CAMphase=np.float64(-1),  # seconds
                          CAMstate=np.int8([-1] * n_piezos),  # =1 if LEDs on for bub t0
                          nPMThits_fastdaq=np.int32([[-1, -1, -1]] * n_piezos),
                          nVetohits_fastdaq=np.int32([[-1, -1, -1]] * n_piezos),
                          t_nearestPMThit=np.float64([0] * n_piezos) + np.nan,
                          t_nearestVetohit=np.float64([0] * n_piezos) + np.nan,
                          t_lastPMThit=np.float64(np.nan),
                          PMTmatch_ix=np.int_([-1] * n_piezos),  # int64 on linux
                          PMTmatch_t0=np.float64([-1] * n_piezos),
                          PMTmatch_lag=np.float64([-1] * n_piezos),
                          PMTmatch_iPMThit=np.int32([-1] * n_piezos),
                          PMTmatch_baseline=np.float64([-1] * n_piezos),
                          PMTmatch_baserms=np.float64([-1] * n_piezos),
                          PMTmatch_area=np.nan +
                                        np.zeros([n_piezos, n_windows], dtype=np.float64),
                          PMTmatch_area_nobs=np.nan +
                                             np.zeros([n_piezos, n_windows], dtype=np.float64),
                          PMTmatch_min=np.nan +
                                       np.zeros([n_piezos, n_windows], dtype=np.float64),
                          PMTmatch_max=np.nan +
                                       np.zeros([n_piezos, n_windows], dtype=np.float64),
                          PMTmatch_pulse_area=np.float64([-1] * n_piezos),
                          PMTmatch_pulse_height=np.float64([-1] * n_piezos),
                          PMTmatch_pulse_tstart=np.float64([-1] * n_piezos),
                          PMTmatch_pulse_tend=np.float64([-1] * n_piezos),
                          PMTmatch_pulse_tpeak=np.float64([-1] * n_piezos),
                          PMTmatch_pulse_t10=np.float64([-1] * n_piezos),
                          PMTmatch_pulse_t90=np.float64([-1] * n_piezos),
                          PMTmatch_coinc=np.int16([-1] * n_piezos),
                          PMTmatch_sat=np.int8([-1] * n_piezos),
                          PMTmatch_nphe=np.float64([-1] * n_piezos),
                          PMTmatch_npeaks=np.int32([-1] * n_piezos),
                          PMTmatch_nclusters=np.int32([-1] * n_piezos),
                          PMTmatch_maxpeak=np.float64([-1] * n_piezos)
                          )

    out = default_output
    try:

        if 'PMT_area' in pmtpa:
            n_windows = pmtpa['PMT_area'].shape[1]
        else:
            n_windows = 2

        n_piezos = aa['bubble_t0'].size
        if len(aa['bubble_t0'].shape) == 0:  # correct numpy's scalar bug
            aa['bubble_t0'] = np.zeros(1, dtype=aa['bubble_t0'].dtype) + aa['bubble_t0']

        if ev['fastDAQ']['loaded']:

            trig_ix = np.int32(np.around(ev['fastDAQ']['caldata']
                                       ['pretrigger_samples'][0]))
            dt = ev['fastDAQ']['caldata']['dt'][0]
            bubble_ix = np.int32(0)
            bubble_t0 = np.float64(0)
            if 'bubble_t0' in aa: # two piezos
                bubble_t0 = aa['bubble_t0']
                bubble_ix = np.int32(np.around((aa['bubble_t0'] / dt) + trig_ix))

            if 'CAMgate' in ev['fastDAQ']:
                led_on = ev['fastDAQ']['CAMgate'] > -0.5
                led_switch = np.diff(np.int8(led_on))
                led_switch_on = led_switch == 1
                led_switch_off = led_switch == -1
                led_switch_on_ix = np.nonzero(led_switch_on)[0] + 1
                led_switch_off_ix = np.nonzero(led_switch_off)[0]

                if (led_switch_on_ix.shape[0] > 1) and\
                        (led_switch_on_ix[0] < trig_ix) and\
                        (led_switch_on_ix[-1] >= trig_ix):

                    last_pretrig_led_on_ix =\
                        led_switch_on_ix[led_switch_on_ix < trig_ix][-1]
                    last_pretrig_led_off_ix =\
                        led_switch_off_ix[led_switch_off_ix >
                                          last_pretrig_led_on_ix][0]
                    first_posttrig_led_on_ix =\
                        led_switch_on_ix[led_switch_on_ix >= trig_ix][0]

                    out['CAMperiod'] = (first_posttrig_led_on_ix -
                                        last_pretrig_led_on_ix) * dt
                    out['CAMduty'] = (last_pretrig_led_off_ix -
                                      last_pretrig_led_on_ix) * dt
                    out['CAMphase'] = (trig_ix - last_pretrig_led_on_ix) * dt
                else:
                    pass

                for i_piezo in np.arange(n_piezos):
                    if (bubble_ix[i_piezo] >= 0) and (bubble_ix[i_piezo] < led_on.shape[0]): # two piezos
                        out['CAMstate'][i_piezo] = np.int8(led_on[bubble_ix[i_piezo]])

            if 'PMTtrig' in ev['fastDAQ']:
                PMTtrig_on = ev['fastDAQ']['PMTtrig'] > 1.5
                PMTtrig_switch_on = PMTtrig_on[1:] * (~PMTtrig_on[:-1])
                PMTtrig_times = ev['fastDAQ']['time'][1:][PMTtrig_switch_on]

                for i_piezo in np.arange(n_piezos):
                    out['nPMThits_fastdaq'][i_piezo,0] =\
                        np.int32(np.sum(PMTtrig_switch_on[:bubble_ix[i_piezo]]))
                    out['nPMThits_fastdaq'][i_piezo,1] =\
                        np.int32(np.sum(PMTtrig_switch_on[bubble_ix[i_piezo]:trig_ix]))
                    out['nPMThits_fastdaq'][i_piezo,2] =\
                        np.int32(np.sum(PMTtrig_switch_on[trig_ix:]))
                if PMTtrig_times.shape[0] > 0:
                    out['t_lastPMThit'] = PMTtrig_times[-1]
                    out['t_nearestPMThit'] = np.zeros(bubble_ix.shape, dtype=np.float64)
                    for i_piezo in np.arange(n_piezos):
                        out['t_nearestPMThit'][i_piezo] =\
                            PMTtrig_times[np.argmin(np.abs(PMTtrig_times -
                                                           bubble_t0[i_piezo] + 4.0e-5))]

            if 'VetoCoinc' in ev['fastDAQ']:
                VetoCoinc_on = ev['fastDAQ']['VetoCoinc'] > 1.5
                VetoCoinc_switch_on = VetoCoinc_on[1:] * (~VetoCoinc_on[:-1])
                VetoCoinc_times = ev['fastDAQ']['time'][1:][VetoCoinc_switch_on]

                for i_piezo in np.arange(n_piezos):
                    out['nVetohits_fastdaq'][i_piezo,0] =\
                        np.int32(np.sum(VetoCoinc_switch_on[:bubble_ix[i_piezo]]))
                    out['nVetohits_fastdaq'][i_piezo,1] =\
                        np.int32(np.sum(VetoCoinc_switch_on[bubble_ix[i_piezo]:trig_ix]))
                    out['nVetohits_fastdaq'][i_piezo,2] =\
                        np.int32(np.sum(VetoCoinc_switch_on[trig_ix:]))
                if VetoCoinc_times.shape[0] > 0:
                    for i_piezo in np.arange(n_piezos):
                        out['t_nearestVetohit'][i_piezo] =\
                            VetoCoinc_times[np.argmin(np.abs(VetoCoinc_times -
                                                             bubble_t0[i_piezo] + 4.0e-5))]

        if ('bubble_t0' in aa) and\
                ('PMT_t0_fastdaq' in pmtpa):
            first_in_cluster = np.ones(pmtpa['PMT_t0_fastdaq'].shape[0],
                                       dtype=bool)
            first_in_cluster[1:] = np.diff(pmtpa['PMT_t0_fastdaq']) > 2.0e-5
            for i_piezo in np.arange(n_piezos):
                if (aa['bubble_t0'][i_piezo] > -1) and\
                    (aa['bubble_t0'][i_piezo] >= pmtpa['PMT_t0_fastdaq'][0]) and\
                    (aa['bubble_t0'][i_piezo] <= (pmtpa['PMT_t0_fastdaq'][-1] + 1e-4)):

                    time_from_bubble = pmtpa['PMT_t0_fastdaq'] - aa['bubble_t0'][i_piezo] + 4.0e-5
                    time_from_bubble[~first_in_cluster] = np.inf
                    out['PMTmatch_ix'][i_piezo] = np.argmin(np.abs(time_from_bubble))
                    out['PMTmatch_t0'][i_piezo] = pmtpa['PMT_t0_fastdaq'][out['PMTmatch_ix'][i_piezo]]
                    out['PMTmatch_lag'][i_piezo] = time_from_bubble[out['PMTmatch_ix'][i_piezo]]
                    out['PMTmatch_iPMThit'][i_piezo] = pmtpa['iPMThit'][out['PMTmatch_ix'][i_piezo]]
                    #ipdb.set_trace()
                    for k in out.keys():
                        if (k[:8] == 'PMTmatch') and\
                                (('PMT' + k[8:]) in pmtpa):
                            out[k][i_piezo] = pmtpa['PMT' + k[8:]][out['PMTmatch_ix'][i_piezo]]
                        elif (k[:8] == 'PMTmatch') and\
                                (('pmt' + k[8:]) in pmtpa):
                            out[k][i_piezo] = pmtpa['pmt' + k[8:]][out['PMTmatch_ix'][i_piezo]]

        return out
    except Exception as e:
        return default_output


if __name__ == "__main__":
    from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE
    from SBCcode.DataHandling.ReadBinary import ReadBlock as RB


    runid = "20170928_5"
    ev = GE("/bluearc/storage/SBC-17-data/{runid}/".format(runid=runid), 22,  "fastDAQ")
    ac = RB("/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output/{runid}/AcousticAnalysis_{runid}.bin".format(runid=runid))
    pmt =  RB("/pnfs/coupp/persistent/grid_output/SBC-17/output/{runid}/PMTpulseAnalysis_{runid}.bin".format(runid=runid))
    timing = TimingAnalysis(ev, aa=ac, pmtpa=pmt)



    pass