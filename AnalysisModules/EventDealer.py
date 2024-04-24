import os
import re
import shutil
import time

import numpy as np
import copy
import numpy.matlib

from .AnalyzeDytran import dytranAnalysis as da
from .EventAnalysis import EventAnalysis as eva
# from SBCcode.AnalysisModules.ImageAnalysis import BubbleFinder
from .AcousticT0 import AcousticAnalysis as aa
from .PMTComprehensiveModule import PMTcm as pmtpa
from .PMTfastDAQalignment import PMTandFastDAQalignment as pmtfda
from .PTData import main as ptd
from .TimingAnalysis import TimingAnalysis as ta
from .ExposureAnalysis import ExposureAnalysis as expa
from ..DataHandling.GetSBCEvent import GetEvent as get_event
from ..DataHandling.WriteBinary import WriteBinaryNtupleFile as wb
from ..DataHandling.ReadBinary import ReadBlock as rb

# For Acoustic t0 test
from ..UserCode.John.NewT0 import calculate_t0 as calculate_t0

def BuildEventList(rundir, first_event=0, last_event=-1):
    # Inputs:
    #   rundir: Directory for the run
    #   first_event: Index of first event
    #   last_event: Index of last_event
    # Outputs: A sorted list of events from rundir
    # Possible Issues: If we ever move to a different naming scheme, this needs to be re-written.
    eventdirlist = os.listdir(rundir)
    eventdirlist = filter(lambda fn: (not re.search('^\d+$', fn) is None) and
                                     os.path.isdir(os.path.join(rundir, fn)),
                          eventdirlist)
    eventdirlist = filter(lambda fn: os.path.exists(os.path.join(rundir,
                                                                 *[fn, 'Event.txt'])), eventdirlist)
    eventlist = np.intp(list(eventdirlist))
    eventlist = eventlist[eventlist >= first_event]
    if last_event >= 0:
        eventlist = eventlist[eventlist <= last_event]
    return np.sort(eventlist)


def ProcessFromReconfile(reconfile, datadir='/exp/e961/data/SBC-17-data',
                         dataset='SBC-2017', recondir='.', process_list=None):
    # Inputs:
    #   reconfile: Location of file that will be used for event list
    #   dataset: Indicator used for filtering which analyses to run
    #   recondir: Location of recon data/where we want to output our binary files
    #   process_list: List of analyses modules to run. example: ["acoustic", "event", ""]
    # Outputs: Nothing. Saves binary files to recondir.
    if process_list is None:
        process_list = []  # This is needed since lists are mutable objects. If you have a default argument
                           # as a mutable object, then the default argument can *change* across multiple
                           # function calls since the argument is created ONCE when the function is defined.
    evlist_source = rb(reconfile)

    if not os.path.isdir(recondir):
        os.mkdir(recondir)

    if dataset == 'SBC-2017':
        loadlist = []
        for process in process_list:
            if process.lower().strip() == 'event':
                loadlist.append('event')
            elif process.lower().strip() == 'images':
                loadlist.append('images')
            elif process.lower().strip() == 'exposure':
                loadlist.append('slowDAQ')
                loadlist.append('event')
            elif process.lower().strip() == 'history':
                loadlist.append('slowDAQ')
            elif any(process.lower().strip() == x for x in ['dytran', 'acoustic']):
                loadlist.append('fastDAQ')
            # elif process.lower().strip() == 'pmtpf':
            #     loadlist.append('PMTtraces')
            elif any(process.lower().strip() == x for x in ['pmtfda', 'pmt', 'timing']):
                loadlist.append('fastDAQ')
                loadlist.append('PMTtraces')
        loadlist = list(set(loadlist))
    else:
        loadlist = ['~']

    exposure_out = []

    exposure_default = expa(None)

    process_list = [p.lower().strip() for p in process_list]
    eventlist = np.concatenate((evlist_source['runid'], evlist_source['ev'][:,None]), axis=1)

    for ev in eventlist:
        runname = str(ev[0])+'_'+str(ev[1])
        rundir = os.path.join(datadir, runname)
        t0 = time.time()
        print('Starting event ' + runname + '/' + str(ev[2]))
        thisevent = get_event(rundir, ev[2], *loadlist)
        print('Time to load event:  '.rjust(35) +
              str(time.time() - t0) + ' seconds')

        if "exposure" in process_list:
            # History analysis
            t1 = time.time()
            if dataset == 'SBC-2017':
                try:
                    exposure_out.append(expa(thisevent))
                except:
                    exposure_out.append(copy.deepcopy(exposure_default))
                exposure_out[-1]['runid'] = ev[:2]
                exposure_out[-1]['ev'] = ev[2:]
            et = time.time() - t1
            print('Exposure analysis:  '.rjust(35) + str(et) + ' seconds')

        print('*** Full event analysis ***  '.rjust(35) +
              str(time.time() - t0) + ' seconds')

    if "exposure" in process_list:
        wb(os.path.join(recondir,
                        'ExposureAnalysis_all.bin'), exposure_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    return


def ProcessSingleRun2(rundir, dataset='SBC-2017', recondir='.', process_list=None):
    # Inputs:
    #   rundir: Location of raw data
    #   dataset: Indicator used for filtering which analyses to run
    #   recondir: Location of recon data/where we want to output our binary files
    #   process_list: List of analyses modules to run. example: ["acoustic", "event", ""]
    # Outputs: Nothing. Saves binary files to recondir.
    if process_list is None:
        process_list = []  # This is needed since lists are mutable objects. If you have a default argument
                           # as a mutable object, then the default argument can *change* across multiple
                           # function calls since the argument is created ONCE when the function is defined.
    runname = os.path.basename(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)
    #run_recondir = os.path.join(recondir, runname)
    run_recondir = recondir

    if not os.path.isdir(run_recondir):
        os.mkdir(run_recondir)

    if dataset == 'SBC-2017':
        loadlist = []
        for process in process_list:
            if process.lower().strip() == 'event':
                loadlist.append('event')
            elif process.lower().strip() == 'images':
                loadlist.append('images')
            elif process.lower().strip() == 'history':
                loadlist.append('slowDAQ')
            elif any(process.lower().strip() == x for x in ['dytran', 'acoustic']):
                loadlist.append('fastDAQ')
            # elif process.lower().strip() == 'pmtpf':
            #     loadlist.append('PMTtraces')
            elif any(process.lower().strip() == x for x in ['pmtfda', 'pmt', 'timing']):
                loadlist.append('fastDAQ')
                loadlist.append('PMTtraces')
        loadlist = list(set(loadlist))
    else:
        loadlist = ['~']

    event_out = []
    pmtfda_out = []
    pmt_out = []
    dytran_out = []
    acoustic_out = []
    history_out = []
    timing_out = []
    image_out = []

    event_default = eva(None)
    pmtfda_default = pmtfda(None)
    pmt_default = pmtpa(None)
    dytran_default = da(None)
    acoustic_default = aa(None, None)
    history_default = ptd(None)
    timing_default = ta(None, None, None)
    # image_default = BubbleFinder(None, None, None ,None, None, None)

    process_list = [p.lower().strip() for p in process_list]
    print("Starting run " + rundir)
    eventlist = BuildEventList(rundir)

    for ev in eventlist:
        t0 = time.time()
        print('Starting event ' + runname + '/' + str(ev))
        npev = np.array([ev], dtype=np.int32)
        thisevent = get_event(rundir, ev, *loadlist)
        print('Time to load event:  '.rjust(35) +
              str(time.time() - t0) + ' seconds')

        # if thisevent["PMTtraces"]["loaded"]:
        #     if thisevent['PMTtraces']['traces'].nbytes > 500 * 1024 ** 2:  # 500MB, skip
        #         thisevent['PMTtraces']['traces'] = np.array([])
        #         thisevent['PMTtraces']['loaded'] = False

        if "event" in process_list:
            # zeroth order of business:  copy event data
            t1 = time.time()
            if dataset == 'SBC-2017':
                try:
                    event_out.append(eva(thisevent))
                except:
                    event_out.append(copy.deepcopy(event_default))
                event_out[-1]['runid'] = runid
                event_out[-1]['ev'] = npev
            et = time.time() - t1
            print('Event analysis:  '.rjust(35) + str(et) + ' seconds')

        if "pmtfda" in process_list:
            # PMTfastDAQalignment
            t1 = time.time()
            if dataset == 'SBC-2017':
                try:
                    pmtfda_out.append(pmtfda(thisevent))
                except:
                    pmtfda_out.append(copy.deepcopy(pmtfda_default))
                pmtfda_out[-1]['runid'] = runid
                pmtfda_out[-1]['ev'] = npev
            et = time.time() - t1
            print('PMT - fastDAQ alignment analysis:  '.rjust(35) + str(et) + ' seconds')

        if "pmt" in process_list:
            t1 = time.time()
            # PMTtraces after pmtfda
            if dataset == 'SBC-2017':
                try:
                    pmt_out.append(pmtpa(thisevent, pmtfda=pmtfda_out[-1]))
                except:
                    pmt_out.append(copy.deepcopy(pmt_default))
                pmt_out[-1]['runid'] = runid + np.zeros((pmt_out[-1][list(pmt_out[-1].keys())[0]].shape[0], 2),
                                                        dtype=np.int32)
                pmt_out[-1]['ev'] = npev + np.zeros(pmt_out[-1][list(pmt_out[-1].keys())[0]].shape[0],
                                                    dtype=np.int32)
            et = time.time() - t1
            print('PMT pulse analysis:  '.rjust(35) + str(et) + ' seconds')

        if "acoustic" in process_list:
            # Acoustic analysis
            t1 = time.time()
            if dataset == 'SBC-2017':
                tau_peak = 0.0025884277467056165  # <-- This comes from TauResultAnalysis.py (J.G.)
                tau_average = 0.0038163479219674467  # <-- This also ^^
                lower_f = 20000
                upper_f = 40000
                piezo_fit_type = 0
                try:
                    acoustic_out.append(aa(ev=thisevent, tau=tau_average, piezo_fit_type=piezo_fit_type,
                                           corr_lowerf=lower_f, corr_upperf=upper_f))
                except:
                    acoustic_out.append(copy.deepcopy(acoustic_default))
                acoustic_out[-1]['runid'] = runid
                acoustic_out[-1]['ev'] = npev
            et = time.time() - t1
            print('Acoustic analysis:  '.rjust(35) + str(et) + ' seconds')

        if "dytran" in process_list:
            # Dytran analysis
            t1 = time.time()
            if dataset == 'SBC-2017':
                try:
                    dytran_out.append(da(thisevent))
                except:
                    dytran_out.append(copy.deepcopy(dytran_default))
                dytran_out[-1]['runid'] = runid
                dytran_out[-1]['ev'] = npev
            et = time.time() - t1
            print('Dytran analysis:  '.rjust(35) + str(et) + ' seconds')

        if "history" in process_list:
            # History analysis
            t1 = time.time()
            if dataset == 'SBC-2017':
                try:
                    history_out.append(ptd(thisevent,
                                           edge=np.arange(13.5, 101.5, 1,
                                                          dtype=np.float64), targetPT='PT6'))
                except:
                    history_out.append(copy.deepcopy(history_default))
                history_out[-1]['runid'] = runid
                history_out[-1]['ev'] = npev
            et = time.time() - t1
            print('History analysis:  '.rjust(35) + str(et) + ' seconds')

        if "timing" in process_list:
            # Timing analysis after acoustic and pmt
            t1 = time.time()
            if dataset == 'SBC-2017':
                try:
                    timing_out.append(ta(thisevent,
                                         pmt_out[-1],
                                         acoustic_out[-1]))
                except:
                    timing_out.append(copy.deepcopy(timing_default))
                timing_out[-1]['runid'] = runid
                timing_out[-1]['ev'] = npev
            et = time.time() - t1
            print('Timing analysis:  '.rjust(35) + str(et) + ' seconds')

        # if "pmtpf" in process_list:
        #     # pmt pulse finder
        #     t1 = time.time()
        #     if dataset == 'SBC-2017':
        #         base_samples = np.intp(80)
        #         amp_gains = np.array([1, 1], dtype=np.float64)
        #
        #         pmtpf_out.append(pmtpf2(thisevent, base_samples, amp_gains=amp_gains))
        #         pmtpf_out[-1]['runid'] = np.matlib.repmat(runid, pmtpf_out[-1]['iTrace'].shape[0], 1)
        #         pmtpf_out[-1]['ev'] = np.matlib.repmat(npev, pmtpf_out[-1]['iTrace'].shape[0], 1)
        #     et = time.time() - t1
        #     print('PMT - pulse finder analysis:  '.rjust(35) + str(et) + ' seconds')

        if False and ("images" in process_list):
            # image analysis
            t1 = time.time()
            if dataset == 'SBC-2017':
                n_pre_trig = 12
                n_frames = 3
                adc_thresh = 15
                bubble_thresh = 4
                try:
                    image_out.append(BubbleFinder(rundir, ev,
                                                  n_pre_trig, n_frames,
                                                  adc_thresh, bubble_thresh))
                except:
                    image_out.append(copy.deepcopy(image_default))
            et = time.time() - t1
            print('Image analysis:  '.rjust(35) + str(et) + ' seconds')

        print('*** Full event analysis ***  '.rjust(35) +
              str(time.time() - t0) + ' seconds')


    #for process in process_list:
    if "event" in process_list:
        wb(os.path.join(run_recondir,
                        'EventAnalysis_' + runname + '.bin'), event_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    if "pmtfda" in process_list:
        wb(os.path.join(run_recondir,
                        'PMTfastDAQalignment_' + runname + '.bin'), pmtfda_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    if "pmt" in process_list:
        wb(os.path.join(run_recondir,
                        'PMTpulseAnalysis_' + runname + '.bin'), pmt_out,
           rowdef=7, initialkeys=['runid', 'ev'], drop_first_dim=True)

    if "dytran" in process_list:
        wb(os.path.join(run_recondir,
                        'DytranAnalysis_' + runname + '.bin'), dytran_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    if "acoustic" in process_list:
        wb(os.path.join(run_recondir,
                        'AcousticAnalysis_' + runname + '.bin'), acoustic_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    if "history" in process_list:
        wb(os.path.join(run_recondir,
                        'HistoryAnalysis_' + runname + '.bin'), history_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    if "timing" in process_list:
        wb(os.path.join(run_recondir,
                        'TimingAnalysis_' + runname + '.bin'), timing_out,
           rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)

    # if "pmtpf" in process_list:
    #     wb(os.path.join(run_recondir,
    #                     'PMTpulseFinder_' + runname + '.bin'), pmtpf_out,
    #        rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)

    if False and "images" in process_list:
        wb(os.path.join(recondir,
                        'ImageAnalysis_' + runname + '.bin'), image_out,
                        rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)

    # copy handscan result
    # if dataset == 'SBC-2015':
    #     human_gbub_dir = '/coupp/data/home/coupp/HumanGetBub_output_SBC-15'
    #     if os.path.isdir(human_gbub_dir):
    #         human_gbub_file = os.path.join(human_gbub_dir,
    #                                        'HumanGetBub_' +
    #                                        runname + '.bin')
    #         if os.path.exists(human_gbub_file):
    #            shutil.copy(human_gbub_file, run_recondir)
    return



if __name__ == "__main__":
    if False:
        ProcessSingleRun2(rundir="/bluearc/storage/SBC-17-data/20171007_3",
                          recondir="/nashome/j/jgresl/Test/Actuals", # Use your own directory for testing~
                          process_list = ["acoustic"])
    elif True:
        ProcessFromReconfile("/pnfs/coupp/persistent/grid_output/SBC-17/output/HistoryAnalysis_all.bin",
                             process_list=["exposure"],
                             recondir="/exp/e961/app/home/coupp/ProgramTesting")
    elif True:
        ProcessFromReconfile("/pnfs/coupp/persistent/grid_output/SBC-17/output/20171011_10/HistoryAnalysis_20171011_10.bin",
                             process_list=["exposure"],
                             recondir="/exp/e961/app/home/coupp/ProgramTesting")
    pass