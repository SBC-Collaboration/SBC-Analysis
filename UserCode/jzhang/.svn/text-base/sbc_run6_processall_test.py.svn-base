# python sbc_pmttest_processall.py [run_list]
# if run_list is provided, the runs in the list will be processed; otherwise
# the runs in the script will be processed

import numpy as np
import SBCcode as sbc
import os
import sys
import socket
import re
import pdb

hostname = socket.gethostname()
if hostname == 'Pegasus':
    datadir = '/mnt/XENON_DAQ/SBC-17-data'
elif re.search('couppgpvm', hostname):
    datadir = '/bluearc/storage/SBC-17-data'

# recondir = '/bluearc/storage/recon/devel/SBC-15/output'

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
    #runlist = ['20170626_5']
    # runlist = ['20170629_1'] # acoustic test
    #runlist = ['20170624_4']
    # runlist = ['20170802_2']
    # runlist = ['20170706_7']
    # runlist = ['20170704_2'] # test dytran
    runlist = ['20170924_1'] # test memory usage
    # runlist = ['20171007_4'] # test memory usage, 1.2G too large
    # runlist = ['20170929_1'] # 500M PMTtraces
    # runlist = ['20170929_0'] # 1.2 G PMTtraces
    runlist = ['20171011_0'] # test PMTPulseAnalysis change
for runname in runlist:
    runid_str = runname.split('_')
    runid = np.int32(runid_str)
    sbc.psr2(os.path.join(datadir, runname),
                  dataset='SBC-2017',
                  recondir=os.path.join(recondir, runname),
            # process_list=['images'])
            #     process_list=['pmtfda','pmt'])
             # process_list=['event', 'pmtfda', 'pmt', 'dytran', 'acoustic', 'history', 'timing', 'images'])
             # process_list=['images'])
    # process_list=['pmtfda', 'pmt', 'acoustic', 'timing'])
    # process_list=['event', 'pmtfda', 'pmt', 'acoustic', 'timing'])
    process_list=['event', 'pmtfda', 'pmt', 'dytran', 'acoustic', 'history', 'timing', 'images'])
    # process_list=['event', 'pmtfda', 'pmt', 'dytran', 'acoustic', 'history', 'timing', 'pmtpf', 'images'])
    # process_list=['event', 'dytran', 'acoustic', 'history', 'images'])
