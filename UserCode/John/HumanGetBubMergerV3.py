
import os
from copy import deepcopy

import numpy as np

from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
from SBCcode.Tools import SBCtools
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as WB

HumanGetBubCache = {}
counter = {"ct": 0, "mask": set()}

def return_HumanGetBub_data(data_path, runid, ev):
    runidstr = str(runid[0]) + "_" + str(runid[1])
    potential_bubble_file = os.path.join(data_path, "HumanGetBub_"+runidstr+".bin")
    default_out = {"runid": np.array([runid]),
                   "ev": np.array([ev]),
                   "ibubimage": np.zeros(1) - 1,
                   "nbubimage": np.zeros(1) - 1,
                   "cam": np.zeros(1) - 1,
                   "frame": np.zeros(1) - 1,
                   "ipix": np.zeros(1) - 1,
                   "jpix": np.zeros(1) - 1
                   }
    if not os.path.isfile(potential_bubble_file):
        counter["ct"] += 1
        return deepcopy(default_out)
    else:
        out = {}
        if runidstr not in HumanGetBubCache.keys():
            HumanGetBubCache[runidstr] = RB(potential_bubble_file)
        mask = HumanGetBubCache[runidstr]["ev"] == ev
        if np.sum(mask) == 0:
            return deepcopy(default_out)
        counter["mask"].add(np.sum(mask))
        for k,v in HumanGetBubCache[runidstr].items():
            out[k] = v[mask]
        return out



if __name__ == "__main__":
    bubble_path = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/"
    output_path = "/pnfs/coupp/persistent/grid_output/SBC-17/output/HumanGetBub_all.bin"
    timing_path = "/pnfs/coupp/persistent/grid_output/SBC-17/output/TimingAnalysis_all.bin"

    # So we want to open our timing data and for each event, see if we have a bubble file. If we have a bubble file,
    # then open it and add it to the dictionary. If not, append a default value.

    do_merge = True
    if do_merge:
        total_out = {}
        first_file = True
        timing_data = RB(timing_path)
        for runid, ev in zip(timing_data["runid"], timing_data["ev"]):
            runidev_str = str(runid[0]) + "_" + str(runid[1]) + "_" + str(ev)
            if first_file:
                total_out = return_HumanGetBub_data(bubble_path, runid, ev)
                first_file = False
            else:
                total_out = SBCtools.dictionary_append(total_out, return_HumanGetBub_data(bubble_path, runid, ev), make_copy=True)

        c1 = total_out["ibubimage"]<2
        if np.all(total_out["ev"][c1] == timing_data["ev"]):
            print("Merge succeeded. Writing to {}".format(output_path))
            WB(output_path, total_out, rowdef=8)
        else:
            print("Merge failed.")