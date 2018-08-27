# python sbc_pmttest_processall.py [run_list]
# if run_list is provided, the runs in the list will be processed; otherwise
# the runs in the script will be processed

import numpy as np
import SBCcode as sbc
import os
import sys

datadir = '/bluearc/storage/SBC-17-data'
# recondir = '/bluearc/storage/recon/devel/SBC-15/output'
# datadir = '/mnt/XENON_DAQ/SBC-17-data'
recondir = '/bluearc/storage/recon/devel/sbc-17/output'

# ~ runlist = os.listdir(datadir)
# ~ runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None) and
# ~                  os.path.isdir(os.path.join(datadir, fn)),
# ~                  runlist)
# ~ runlist = filter(lambda fn: os.path.exists(os.path.join(datadir,
# ~                 *[fn, 'DAQversion.txt'])), runlist)


if len(sys.argv) > 1:
    #runlist = sys.argv[1:]
    istart = int(sys.argv[1])
    iend = int(sys.argv[2])
else:
    istart = 0
    iend = np.iinfo(np.int32).max

runlist=[]
with open('list_20170705','r') as f:
    runlist = f.readlines()
runlist = [run.strip() for run in runlist]

#for runname in runlist:
for ii in range(istart-1, np.amin([iend, len(runlist)])):
    runname = runlist[ii]
    runid_str = runname.split('_')
    runid = np.int32(runid_str)
    sbc.psr2(os.path.join(datadir, runname),
                  dataset='SBC-2017',
                  recondir=os.path.join(recondir, runname),
                  process_list=['event', 'pmtfda', 'pmt', 'dytran', 'acoustic', 'history', 'timing', 'images'])
