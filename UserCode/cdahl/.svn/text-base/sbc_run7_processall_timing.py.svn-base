import numpy as np
import SBCcode as sbc
import os
import re
import time
import copy
import sys
import gc
from SBCcode.AnalysisModules.EventDealer import BuildEventList as bevl
from SBCcode.AnalysisModules.TimingAnalysis import TimingAnalysis as ta
from SBCcode.DataHandling.GetSBCEvent import GetEvent as get_event
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as wb

datadir = '/bluearc/storage/SBC-17-data'
recondir_pmt = '/pnfs/coupp/persistent/grid_output/SBC-17/output'
recondir_aa = '/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output'
recondir_xyz = '/pnfs/coupp/persistent/grid_output/SBC-17/output'
recondir_output = '/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output'

runlist = os.listdir(recondir_pmt)
runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None) and
                 os.path.isdir(os.path.join(datadir, fn)),
                 runlist)
runlist = filter(lambda fn: os.path.exists(os.path.join(recondir_pmt,
                 *[fn, 'PMTpulseAnalysis_' +
                   fn + '.bin'])), runlist)
runlist = filter(lambda fn: os.path.exists(os.path.join(recondir_aa,
                 *[fn, 'AcousticAnalysis_' +
                   fn + '.bin'])), runlist)
runlist = filter(lambda fn: not os.path.exists(os.path.join(recondir_aa,
                 *[fn, 'TimingAnalysis_' +
                   fn + '.bin'])), runlist)
xyz = sbc.read_bin(os.path.join(recondir_xyz, 'SimpleXYZ_all.bin'))
runlist = filter(lambda fn: np.any(np.all(xyz['runid']==np.int32(fn.split('_')),axis=1)), runlist)

timing_default = ta(None, None, None)

for runname in runlist:
    rundir = os.path.join(datadir, runname)
    evlist = bevl(rundir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    timing_out = []
    pmtpa = None
    gc.collect()

    try:
        pmtpa = sbc.read_bin(os.path.join(recondir_pmt, *[runname, 'PMTpulseAnalysis_' + runname + '.bin']), 5000)
    except:
        print('could not load pmtpa')
        sys.stdout.flush()
        continue

    aa = sbc.read_bin(os.path.join(recondir_aa, *[runname, 'AcousticAnalysis_' + runname + '.bin']))

    for ev in evlist:
        t0 = time.time()
        print('Starting event ' + runname + '/' + str(ev))
        npev = np.array([ev], dtype=np.int32)
        try:
            thisevent = get_event(rundir, ev, 'fastDAQ')
        except:
            timing_out.append(copy.deepcopy(timing_default))
            timing_out[-1]['runid'] = runid
            timing_out[-1]['ev'] = npev
            print('Failed to load event, skpping')
            sys.stdout.flush()
            continue
        print('Time to load event:  '.rjust(35) +
              str(time.time() - t0) + ' seconds')
        sys.stdout.flush()

        xyzcut = np.all(xyz['runid'] == runid, axis=1) * (xyz['ev'] == npev)
        if np.any(xyzcut):
            thisxyz = dict(bubZ=xyz['bubZ'][xyzcut][0])
        else:
            thisxyz = []
            print('Whoa, no xyz for this event')

        aacut = np.all(aa['runid'] == runid, axis=1) * (aa['ev'] == npev)
        thisaa = dict(bubble_t0=aa['bubble_t0'][aacut][0])

        pmtpacut = np.all(pmtpa['runid'] == runid, axis=1) * (pmtpa['ev'] == npev)

        t1 = time.time()
        timing_out.append(ta(thisevent,
                             pmtpa,
                             thisaa,
                             thisxyz,
                             pmtpacut))

        timing_out[-1]['runid'] = runid
        timing_out[-1]['ev'] = npev
        et = time.time() - t1
        print('Timing analysis:  '.rjust(35) + str(et) + ' seconds')
        sys.stdout.flush()

    wb(os.path.join(recondir_aa,
                    *[runname, 'TimingAnalysis_' + runname + '.bin']), timing_out,
       rowdef=1, initialkeys=['runid', 'ev'], drop_first_dim=False)
    print('Wrote file TimingAnalysis_' + runname + '.bin')
    sys.stdout.flush()
print('All done!')

