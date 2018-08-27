# python sbc_pmttest_processall.py [run_list]
# if run_list is provided, the runs in the list will be processed; otherwise
# the runs in the script will be processed

import SBCcode.AnalysisModules.ImageAnalysis as ia
import SBCcode.DataHandling.WriteBinary as wb
import numpy as np
# import SBCcode as sbc
import os
import re
import sys

# datadir = '/bluearc/storage/SBC-17-data'
#recondir = '/bluearc/storage/recon/devel/SBC-15/output'
datadir = '/mnt/XENON_DAQ/SBC-17-data'
recondir = '.'

# ~ runlist = os.listdir(datadir)
# ~ runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None) and
# ~                  os.path.isdir(os.path.join(datadir, fn)),
# ~                  runlist)
# ~ runlist = filter(lambda fn: os.path.exists(os.path.join(datadir,
# ~                 *[fn, 'DAQversion.txt'])), runlist)

if len(sys.argv) > 1:
    runlist = sys.argv[1:]
else:
    runlist = ['20170625_0']


for runname in runlist:
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    rundir = os.path.join(datadir,runname)
    eventdirlist = os.listdir(rundir)
    eventdirlist = filter(lambda fn: (not re.search('^\d+$', fn) is None) and
                                     os.path.isdir(os.path.join(rundir, fn)),
                          eventdirlist)
    eventlist = [int(x) for x in list(eventdirlist)]
    eventlist = [21]
    if not os.path.isdir(recondir):
        os.mkdir(recondir)
    bubbleList = []
    for ev in eventlist:
        bubbleList.append(ia.BubbleFinder(os.path.join(datadir,runname), ev,
                     12, 3, 15, 4).bubbles)
    #print(bubbleList)
    wb.WriteBinaryNtupleFile(os.path.join(recondir,'ImageAnalysis_' + runname + '.bin'), bubbleList,
        rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
