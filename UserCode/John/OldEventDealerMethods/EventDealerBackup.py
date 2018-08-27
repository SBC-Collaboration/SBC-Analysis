'''from .AlternateAcousticT0 import T0finder as aa
from .AlternateAcousticT0 import T0finder2 as aa2
from .AnalyzeDytran import dytranAnalysis as da
from .PTData import main as ptd
from .PMTfastDAQalignment import PMTandFastDAQalignment as pmtfda
from .TimingAnalysis import TimingAnalysis as ta
from .EventAnalysis import EventAnalysis as eva
from .PMTphe_counter import PMTphea as pmtphea
from .PMTPulseFinder import PMTPulseFinder as pmtpf
from .PMTPulseFinder import PMTPulseFinder2 as pmtpf2
from .ImageAnalysis import BubbleFinder
from .AltAlternateAcousticT0 import run as corr_t0_finder





def ProcessSingleRun_pmtfdaonly(rundir, dataset='SBC-2015', recondir='.'):
     The usual ProcessSingleRun script

        This will do the full analysis for a given run, pmtfda only,
        what it does is determined by the chamber argument,
        where it puts it is determined by the recondir argument


    eventlist = BuildEventList(rundir)

    runname = os.path.basename(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    if False:
        run_recondir = os.path.join(recondir, runname)
    else:
        run_recondir = recondir

    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    if dataset == 'SBC-2015':
        loadlist = ['fastDAQ', 'PMTtraces', 'event']
    else:
        loadlist = ['~']

    pmtfda_out = []

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + str(ev))
        npev = np.array([ev], dtype=np.int32)
        thisevent = get_event(rundir, ev, *loadlist)
        print('Time to load event:  ' +
              str(time.time() - t0) + ' seconds')

        # PMTfastDAQalignment only
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            pmtfda_out.append(pmtfda(thisevent))
            pmtfda_out[-1]['runid'] = runid
            pmtfda_out[-1]['ev'] = npev
        et = time.time() - t1
        print('PMT - fastDAQ alignment analysis:  ' + str(et) + ' seconds')

        print('*** Full event analysis ***  ' +
              str(time.time() - t0) + ' seconds')

    if dataset == 'SBC-2015':
        # pdb.set_trace()
        wb(os.path.join(run_recondir,
                        'PMTfastDAQalignment_' + runname + '.bin'), pmtfda_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    return




def ProcessSingleRun_phecountingonly(rundir, dataset='SBC-2015', recondir='.'):
    The usual ProcessSingleRun script

        This will do the full analysis for a given run, phe counting only,
        what it does is determined by the chamber argument,
        where it puts it is determined by the recondir argument


    eventlist = BuildEventList(rundir)

    runname = os.path.basename(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    if False:
        run_recondir = os.path.join(recondir, runname)
    else:
        run_recondir = recondir

    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    if dataset == 'SBC-2015' or dataset == 'SBC-2017':
        loadlist = ['PMTtraces']
    else:
        loadlist = ['~']

    phe_out = []

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + os.path.join(rundir, str(ev)))
        npev = np.array([ev], dtype=np.int32)
        thisevent = get_event(rundir, ev, *loadlist)
        print('Time to load event:  ' +
              str(time.time() - t0) + ' seconds')

        # PHE Analysis only
        t1 = time.time()
        if dataset == 'SBC-2015' or dataset =='SBC-2017' and True:
            phe_out.append(pmtphea(thisevent))
            phe_out[-1]['runid'] = runid +\
                np.zeros((phe_out[-1][list(phe_out[-1].keys()
                                           )[0]].shape[0], 2),
                         dtype=np.int32)
            phe_out[-1]['ev'] = npev +\
                np.zeros(phe_out[-1][list(phe_out[-1].keys())[0]].shape[0],
                         dtype=np.int32)
        et = time.time() - t1
        print('PHE analysis:  ' + str(et) + ' seconds')

        print('*** Full event analysis ***  ' +
              str(time.time() - t0) + ' seconds')

    if dataset == 'SBC-2015' or dataset == 'SBC-2017':
        # pdb.set_trace()
        wb(os.path.join(run_recondir,
                        'PMTpheAnalysis_' + runname + '.bin'), phe_out,
           rowdef=7, initialkeys=['runid', 'ev'], drop_first_dim=True)

    return



def ProcessSingleRun_historyonly(rundir, dataset='SBC-2015', recondir='.'):
    The usual ProcessSingleRun script

        This will do the full analysis for a given run, pmtfda only,
        what it does is determined by the chamber argument,
        where it puts it is determined by the recondir argument


    eventlist = BuildEventList(rundir)

    runname = os.path.basename(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    if False:
        run_recondir = os.path.join(recondir, runname)
    else:
        run_recondir = recondir

    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    if dataset == 'SBC-2015':
        loadlist = ['slowDAQ']
    else:
        loadlist = ['~']

    history_out = []

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + str(ev))
        npev = np.array([ev], dtype=np.int32)
        thisevent = get_event(rundir, ev, *loadlist)
        print('Time to load event:  ' +
              str(time.time() - t0) + ' seconds')

        # History Analysis only
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            history_out.append(ptd(thisevent,
                                   edge=np.arange(19.5, 101.5, 1,
                                                  dtype=np.float64)))
            history_out[-1]['runid'] = runid
            history_out[-1]['ev'] = npev
        et = time.time() - t1
        print('History analysis:  ' + str(et) + ' seconds')

        print('*** Full event analysis ***  ' +
              str(time.time() - t0) + ' seconds')

    if dataset == 'SBC-2015':
        # pdb.set_trace()
        wb(os.path.join(run_recondir,
                        'HistoryAnalysis_' + runname + '.bin'), history_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    return


def ProcessSingleRun(rundir, dataset='SBC-2015', recondir='.'):
    The usual ProcessSingleRun script

        This will do the full analysis for a given run,
        what it does is determined by the chamber argument,
        where it puts it is determined by the recondir argument


    eventlist = BuildEventList(rundir)

    runname = os.path.basename(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    if False:
        run_recondir = os.path.join(recondir, runname)
    else:
        run_recondir = recondir

    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    if dataset == 'SBC-2015':
        loadlist = ['fastDAQ', 'slowDAQ', 'PMTtraces', 'event']
    else:
        loadlist = ['~']

    event_out = []
    pmtfda_out = []
    pmt_out = []
    dytran_out = []
    acoustic_out = []
    history_out = []
    timing_out = []

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + str(ev))
        npev = np.array([ev], dtype=np.int32)
        thisevent = get_event(rundir, ev, *loadlist)
        print('Time to load event:  ' +
              str(time.time() - t0) + ' seconds')

        # zeroth order of business:  copy event data
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            event_out.append(eva(thisevent))
            event_out[-1]['runid'] = runid
            event_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Event analysis:  ' + str(et) + ' seconds')

        # PMTfastDAQalignment first
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            pmtfda_out.append(pmtfda(thisevent))
            pmtfda_out[-1]['runid'] = runid
            pmtfda_out[-1]['ev'] = npev
        et = time.time() - t1
        print('PMT - fastDAQ alignment analysis:  ' + str(et) + ' seconds')

        # PMTtraces after pmtfda
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            pmt_out.append(pmtpa(thisevent,
                           pmtfda=pmtfda_out[-1]))
            pmt_out[-1]['runid'] = runid +\
                np.zeros((pmt_out[-1][list(pmt_out[-1].keys()
                                           )[0]].shape[0], 2),
                         dtype=np.int32)
            pmt_out[-1]['ev'] = npev +\
                np.zeros(pmt_out[-1][list(pmt_out[-1].keys())[0]].shape[0],
                         dtype=np.int32)
        et = time.time() - t1
        print('PMT pulse analysis:  ' + str(et) + ' seconds')

        # acoustic analysis
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            acoustic_out.append(aa(thisevent))
            acoustic_out[-1]['runid'] = runid
            acoustic_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Acoustic analysis:  ' + str(et) + ' seconds')

        # dytran analysis
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            dytran_out.append(da(thisevent))
            dytran_out[-1]['runid'] = runid
            dytran_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Dytran analysis:  ' + str(et) + ' seconds')

        # timing analysis after acoustic and pmt
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            timing_out.append(ta(thisevent,
                                 pmt_out[-1],
                                 acoustic_out[-1]))
            timing_out[-1]['runid'] = runid
            timing_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Timing analysis:  ' + str(et) + ' seconds')

        # history analysis
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            history_out.append(ptd(thisevent,
                                   edge=np.arange(19.5, 101.5, 1,
                                                  dtype=np.float64)))
            history_out[-1]['runid'] = runid
            history_out[-1]['ev'] = npev
        et = time.time() - t1
        print('History analysis:  ' + str(et) + ' seconds')

        print('*** Full event analysis ***  ' +
              str(time.time() - t0) + ' seconds')

    if dataset == 'SBC-2015':
        # pdb.set_trace()
        wb(os.path.join(run_recondir,
                        'EventAnalysis_' + runname + '.bin'), event_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        wb(os.path.join(run_recondir,
                        'PMTfastDAQalignment_' + runname + '.bin'), pmtfda_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        wb(os.path.join(run_recondir,
                        'PMTpulseAnalysis_' + runname + '.bin'), pmt_out,
           rowdef=7, initialkeys=['runid', 'ev'], drop_first_dim=True)
        wb(os.path.join(run_recondir,
                        'AcousticAnalysis_' + runname + '.bin'), acoustic_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        wb(os.path.join(run_recondir,
                        'TimingAnalysis_' + runname + '.bin'), timing_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        wb(os.path.join(run_recondir,
                        'DytranAnalysis_' + runname + '.bin'), dytran_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        wb(os.path.join(run_recondir,
                        'HistoryAnalysis_' + runname + '.bin'), history_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

        if dataset == 'SBC-2015' and True:
            human_gbub_dir = '/coupp/data/home/coupp/HumanGetBub_output_SBC-15'
            if os.path.isdir(human_gbub_dir):
                human_gbub_file = os.path.join(human_gbub_dir,
                                               'HumanGetBub_' +
                                               runname + '.bin')
                if os.path.exists(human_gbub_file):
                    shutil.copy(human_gbub_file, run_recondir)

    return



def ProcessEventList(eventlist_dict, datadir,
                     dataset='SBC-2015', recondir='.'):
    The usual ProcessEventList script

        This will do the full analysis for a list of events,
        given in the eventlist_dict dictionary of runid and ev,
        same as in the output binaries.
        What it does is determined by the chamber argument,
        where it puts it is determined by the recondir argument


    eventlist = eventlist_dict['ev']

    if dataset == 'SBC-2015':
        loadlist = ['slowDAQ']
    else:
        loadlist = ['~']

    pmt_out = []
    dytran_out = []
    acoustic_out = []
    history_out = []
    timing_out = []

    for i_ev in range(eventlist.size):
        ev = eventlist[i_ev]
        runid = eventlist_dict['runid'][i_ev, :]
        rundir = os.path.join(datadir,
                              str(runid[0]) + '_' + str(runid[1]))

        t0 = time.time()
        print('Starting event ' + str(ev))
        npev = np.array([ev], dtype=np.int32)
        thisevent = get_event(rundir, ev, *loadlist)
        print('Time to load event:  ' +
              str(time.time() - t0) + ' seconds')

        # PMTtraces first
        t1 = time.time()
        if dataset == 'SBC-2015' and False:
            pmt_out.append(pmtpa(thisevent))
            pmt_out[-1]['runid'] = runid +\
                np.zeros((pmt_out[-1][list(pmt_out[-1].keys()
                                           )[0]].shape[0], 2),
                         dtype=np.int32)
            pmt_out[-1]['ev'] = npev +\
                np.zeros(pmt_out[-1][list(pmt_out[-1].keys())[0]].shape[0],
                         dtype=np.int32)
        et = time.time() - t1
        print('PMT pulse analysis:  ' + str(et) + ' seconds')

        # acoustic analysis
        t1 = time.time()
        if dataset == 'SBC-2015' and False:
            acoustic_out.append(aa(thisevent))
            acoustic_out[-1]['runid'] = runid
            acoustic_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Acoustic analysis:  ' + str(et) + ' seconds')

        # dytran analysis
        t1 = time.time()
        if dataset == 'SBC-2015' and False:
            dytran_out.append(da(thisevent))
            dytran_out[-1]['runid'] = runid
            dytran_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Dytran analysis:  ' + str(et) + ' seconds')

        # timing analysis
        t1 = time.time()
        if dataset == 'SBC-2015' and False:
            timing_out.append(pmtfda(thisevent))
            timing_out[-1]['runid'] = runid
            timing_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Timing analysis:  ' + str(et) + ' seconds')

        # history analysis
        t1 = time.time()
        if dataset == 'SBC-2015' and True:
            history_out.append(ptd(thisevent))
            history_out[-1]['runid'] = runid
            history_out[-1]['ev'] = npev
        et = time.time() - t1
        print('History analysis:  ' + str(et) + ' seconds')

        print('*** Full event analysis ***  ' +
              str(time.time() - t0) + ' seconds')

    if dataset == 'SBC-2015':
        # pdb.set_trace()
        # wb(os.path.join(recondir, 'PMTpulseAnalysis.bin'), pmt_out,
        #    rowdef=7, initialkeys=['runid', 'ev'], drop_first_dim=True)
        # wb(os.path.join(recondir, 'AcousticAnalysis.bin'), acoustic_out,
        #    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        # wb(os.path.join(recondir, 'TimingAnalysis.bin'), timing_out,
        #    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        # wb(os.path.join(recondir, 'DytranAnalysis.bin'), dytran_out,
        #    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
        wb(os.path.join(recondir, 'HistoryAnalysis.bin'), history_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    return (acoustic_out, timing_out, dytran_out)



def ProcessSingleRun_pmtpfonly(rundir, dataset='SBC-2017', recondir='.'):
     The usual ProcessSingleRun script

        This will do the full analysis for a given run, pmtpa only,
        what it does is determined by the chamber argument,
        where it puts it is determined by the recondir argument


    eventlist = BuildEventList(rundir)

    runname = os.path.basename(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    if False:
        run_recondir = os.path.join(recondir, runname)
    else:
        run_recondir = recondir

    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    base_samples=np.intp(80)
    amp_gains = np.array([1, 1], dtype=np.float64)


    if dataset == 'SBC-2015':
        loadlist = ['PMTtraces']
    elif dataset == 'SBC-2017':
        loadlist = ['PMTtraces']
        base_samples=np.intp(100)
        amp_gains = np.array([1, 1], dtype=np.float64)
    else:
        loadlist = ['~']

    pmtpa_out = []
    # eventlist = [0]
    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + os.path.join(rundir, str(ev)))
        npev = np.array([ev], dtype=np.int32)
        thisevent = get_event(rundir, ev, *loadlist)
        print('Time to load event:  ' +
              str(time.time() - t0) + ' seconds')

        # PMTfastDAQalignment only
        t1 = time.time()
        if True:
            pmtpa_out.append(pmtpf2(thisevent, base_samples, amp_gains=amp_gains))
            pmtpa_out[-1]['runid'] = np.matlib.repmat(runid, pmtpa_out[-1]['iTrace'].shape[0], 1)
            pmtpa_out[-1]['ev'] = np.matlib.repmat(npev, pmtpa_out[-1]['iTrace'].shape[0], 1)
        et = time.time() - t1
        print('PMT - pulse integration analysis:  ' + str(et) + ' seconds')

        print('*** Full event analysis ***  ' +
              str(time.time() - t0) + ' seconds')

    if True:
        #pdb.set_trace()
        wb(os.path.join(run_recondir,
                        'PMTPulseAnalysis_' + runname + '.bin'), pmtpa_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)

    return


def RunLevel1AnalysisModules(modulename, rundir, eventlist=None,
                             datafilelist=None, outputdir='.',
                             **moduleparams):
    if eventlist is None:
        eventlist = BuildEventList(rundir)
    if datafilelist is None:
        datafilelist = ['~']

    output_dict = [modulename(get_event(rundir, ev, *datafilelist),
                   **moduleparams) for ev in eventlist]
    wb(output_dict)
    return'''