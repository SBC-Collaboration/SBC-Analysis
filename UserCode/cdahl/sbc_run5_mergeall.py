import numpy as np
import SBCcode as sbc
import os
import re
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as wb
import pdb

event_out = []
pmtfda_out = []
pmt_out = []
phe_out = []
dytran_out = []
acoustic_out = []
history_out = []
timing_out = []
gbub_out = []

# recondir = '/bluearc/storage/recon/devel/SBC-15/output'
recondir = '/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output'

runlist = os.listdir(recondir)
runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None) and
                 os.path.isdir(os.path.join(recondir, *[fn, 'TimingAnalysis_' + fn + '.bin'])),
                 runlist)
try:
    for runname in runlist:
        runid_str = runname.split('_')
        runid = np.int32(runid_str)
        if (runid[0] >= 20161031) and\
                (runid[0] <= 20191107):
            # event_out.append(sbc.read_bin(os.path.join(recondir,
            #                                            runname,
            #                                            'EventAnalysis_' +
            #                                            runname + '.bin')))
            # pmtfda_out.append(sbc.read_bin(os.path.join(recondir,
            #                                             runname,
            #                                             'PMTfastDAQalignment_' +
            #                                             runname + '.bin')))
            # pmt_out.append(sbc.read_bin(os.path.join(recondir,
            #                                          runname,
            #                                          'PMTpulseAnalysis_' +
            #                                          runname + '.bin')))
            # phe_out.append(sbc.read_bin(os.path.join(recondir,
            #                                          runname,
            #                                          'PMTpheAnalysis_' +
            #                                          runname + '.bin')))
            # dytran_out.append(sbc.read_bin(os.path.join(recondir,
            #                                             runname,
            #                                             'DytranAnalysis_' +
            #                                             runname + '.bin')))
            # acoustic_out.append(sbc.read_bin(os.path.join(recondir,
            #                                               runname,
            #                                               'AcousticAnalysis_' +
            #                                               runname + '.bin')))
            # history_out.append(sbc.read_bin(os.path.join(recondir,
            #                                              runname,
            #                                              'HistoryAnalysis_' +
            #                                              runname + '.bin')))
            timing_out.append(sbc.read_bin(os.path.join(recondir,
                                                        runname,
                                                        'TimingAnalysis_' +
                                                        runname + '.bin')))
            # gbub_out.append(sbc.read_bin(os.path.join(recondir,
            #                                           runname,
            #                                           'HumanGetBub_' +
            #                                           runname + '.bin')))
except:
    pdb.set_trace()

runname = 'all'
# wb(os.path.join(recondir,
#                 'EventAnalysis_' + runname + '.bin'), event_out,
#    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
# wb(os.path.join(recondir,
#                 'PMTfastDAQalignment_' + runname + '.bin'), pmtfda_out,
#    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
# wb(os.path.join(recondir,
#                 'PMTpulseAnalysis_' + runname + '.bin'), pmt_out,
#    rowdef=7, initialkeys=['runid', 'ev'], drop_first_dim=True)
# wb(os.path.join(recondir,
#                 'PMTpheAnalysis_' + runname + '.bin'), phe_out,
#    rowdef=7, initialkeys=['runid', 'ev'], drop_first_dim=True)
# wb(os.path.join(recondir,
#                 'AcousticAnalysis_' + runname + '.bin'), acoustic_out,
#    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
wb(os.path.join(recondir,
                'TimingAnalysis_' + runname + '.bin'), timing_out,
   rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
# wb(os.path.join(recondir,
#                 'DytranAnalysis_' + runname + '.bin'), dytran_out,
#    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
# wb(os.path.join(recondir,
#                 'HistoryAnalysis_' + runname + '.bin'), history_out,
#    rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=True)
# wb(os.path.join(recondir,
#                 'HumanGetBub_' + runname + '.bin'), gbub_out,
#    rowdef=8, initialkeys=['runid', 'ev'], drop_first_dim=True)
