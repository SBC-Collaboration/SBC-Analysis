import numpy as np
import SBCcode as sbc
import os
import re

datadir = '/bluearc/storage/SBC-15-data'
#recondir = '/bluearc/storage/recon/devel/SBC-15/output'
recondir = '/scratch/coupp/jzhang/sbc'

runlist = os.listdir(datadir)
runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None) and
                 os.path.isdir(os.path.join(datadir, fn)),
                 runlist)
runlist = filter(lambda fn: os.path.exists(os.path.join(datadir,
                 *[fn, 'DAQversion.txt'])), runlist)

runlist=['20161107_23']
for runname in runlist:
    runid_str = runname.split('_')
    runid = np.int32(runid_str)
    if (runid[0] >= 20161031) and\
            (runid[0] <= 20161107) and\
            not(runid[0] == 20161102 and runid[1] == 7) and\
            not(runid[0] == 20161107 and runid[1] == 21):
        try:
            sbc.psr(os.path.join(datadir, runname),
                    recondir=os.path.join(recondir, runname))
        except:
            print('Failed on run ' + runname)
