# python sbc_pmttest_processall.py [run_list]
# if run_list is provided, the runs in the list will be processed; otherwise
# the runs in the script will be processed

import numpy as np
import SBCcode as sbc
import os
import re
import sys
import pdb

#datadir = '/bluearc/storage/SBC-15-data'
#recondir = '/bluearc/storage/recon/devel/SBC-15/output'
datadir = '/mnt/XENON_DAQ/SBC-17-data'
recondir = '.'


#~ runlist = os.listdir(datadir)
#~ runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None) and
#~                  os.path.isdir(os.path.join(datadir, fn)),
#~                  runlist)
#~ runlist = filter(lambda fn: os.path.exists(os.path.join(datadir,
#~                 *[fn, 'DAQversion.txt'])), runlist)


if len(sys.argv) > 1:
    runlist = sys.argv[1:]
else :
    # runlist = ['20170330_2']
    # runlist = ['20170501_0']
    # runlist = ['20170407_2']
    #runlist = ['20170330_10',
    #'20170330_11',
    #'20170406_0',
    #'20170406_2',
    #'20170406_3',
    #'20170407_1',
    #'20170407_2',
    #'20170407_3',
    #'20170407_8',
    #'20170407_9',
    #'20170407_10',
    #'20170407_11',
    #'20170412_2',
    #'20170412_3',
    #'20170412_4',
    #'20170412_5',
    #'20170412_6',
    #'20170412_7',
    #'20170425_4',
    #'20170425_5',
    #'20170501_0',
    #'20170501_1',
    #'20170501_2',
    #'20170501_3',
    #'20170504_0',
    #'20170504_1',
    #'20170504_2',
    #'20170504_3']
    
    runlist = ['20170330_0',
    '20170330_1',
    '20170330_3',
    '20170330_4',
    '20170330_5',
    '20170330_6',
    '20170330_7',
    '20170330_8',
    '20170330_9',
    '20170407_4',
    '20170407_5',
    '20170407_6',
    '20170407_7',
    '20170425_0',
    '20170425_1',
    '20170425_2',
    '20170425_3',
    '20170512_2 ',
    '20170518_0']
    runlist = ['20170512_2']


for runname in runlist:
    runid_str = runname.split('_')
    runid = np.int32(runid_str)
    if True:
        #try:
      sbc.psr_pmtpf(os.path.join(datadir, runname),
                      dataset='SBC-2017',
                      recondir=os.path.join(recondir, runname))
        #except:
        #    print('Failed on run ' + runname)
