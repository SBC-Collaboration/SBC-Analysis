import numpy as np
import SBCcode as sbc
import os
import re
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as wb
import pdb

recondir = '/bluearc/storage/recon/devel/SBC-15/output'

dd = sbc.read_bin(os.path.join(recondir, 'PMTpulseAnalysis_all.bin'))
ev_float = np.float64(dd['runid'][:, 0] - 20161000) +\
    np.float64(dd['runid'][:, 1]) * 1e-3 +\
    np.float64(dd['ev']) * 1e-6

firsthit = dd['iPMThit'] < 2
lasthit = np.ones(firsthit.shape, dtype=np.bool)
lasthit[:-1] = firsthit[1:]

dd_oneline = {k: dd[k][lasthit] for k in dd.keys()}

wb(os.path.join(recondir, 'PMTpulseAnalysis_all_oneline.bin'), [dd_oneline],
    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
