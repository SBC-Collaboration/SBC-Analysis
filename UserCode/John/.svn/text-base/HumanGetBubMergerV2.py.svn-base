
import os
from copy import deepcopy
from collections import defaultdict, OrderedDict

import numpy as np

from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as WB

def sort_bubble_files(arr):
    # Input:
    #   arr: An array of bubble file names as strings. Should look like ["20170623_0", "20170623_5", etc...]
    # Outputs: A natural-ish sorted version that puts the dates in order and the run numbers for each date in order
    # Modified from sort_runs()
    dates_only = []
    runs_only = []
    for bub_id in arr:
        bub_id = bub_id.split(".")[0] # Everything before the .bin
        dates_only.append(bub_id.split("_")[1])
        runs_only.append(bub_id.split("_")[2])
    run_dict = defaultdict(list)
    for date, run in zip(dates_only, runs_only):
        run_dict[date].append(run)
    k = sorted(list(run_dict.keys()), key=int)
    out_list = []
    for date in k:
        run_ids_d = sorted(run_dict[date], key=int)
        for run_id in run_ids_d:
            out_list.append(date+"_"+run_id)
    return out_list

def dictionary_append(d1, *args):
    # Inputs:
    #   d1: Dictionary
    #   args: Any number of dictionaries
    # Outputs: A dictionary with the same keys, but the values are the values from args appended to the values of d1
    # Note: The keys in d1 and args MUST match, and all of the values MUST be lists.
    d1 = deepcopy(d1)
    if len(args) == 0:
        raise TypeError("dictionary_append must be called with at least 2 arguments.")
    for arg in args:
        if type(arg) not in [dict, defaultdict, OrderedDict]:
            raise TypeError("args must be dictionaries!")
    for d2 in args:
        if not set(d1.keys()) == set(d2.keys()):
            raise KeyError("The keys for the two dictionaries must match! Mismatched keys = {}".\
                           format(set(d1.keys()).symmetric_difference(set(d2.keys()))))
    for d2 in args: # This '2nd' for loop because we want to make sure all the keys are the same first.
        for k,v in d2.items():
            if type(v) not in [list, np.ndarray]:
                raise ValueError("The values of the dictionary MUST be list or np.ndarray. Key {} has type(value)={}".\
                             format(k, type(v)))
            try:
                d1[k].extend(v)   # <-- If we have python lists
            except AttributeError:
                # print("DEBUG:", k, d1[k].shape, v.shape)
                d1[k] = np.append(d1[k], v, axis=0) # <-- If we have numpy arrays
    return d1




if __name__ == "__main__":
    bub_data_path = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/"
    #output_path = "/pnfs/coupp/persistent/grid_output/SBC-17/output/"
    timing_path = "/pnfs/coupp/persistent/grid_output/SBC-17/output/TimingAnalysis_all.bin"

    output_path = "/nashome/j/jgresl/BubbleMergeAll.bin"

    default_output = OrderedDict()
    default_output["ibubimage"] = np.zeros(1) - 1.
    default_output["nbubimage"] = np.zeros(1) - 1.
    default_output["cam"] = np.zeros(1) - 1.
    default_output["frame"] = np.zeros(1) - 1.
    default_output["ipix"] = np.zeros(1) - 1.
    default_output["jpix"] = np.zeros(1) - 1.
    do_merge = True
    if do_merge:
        timing_data = RB(timing_path)
        bubble_files = sort_bubble_files([b for b in os.listdir(bub_data_path) if\
                                          os.path.isfile(os.path.join(bub_data_path, b))])
        print(len(bubble_files))
        out = {}
        first_file = True


        prev_runid = ""

        for runid, ev in zip(timing_data["runid"], timing_data["ev"]):
            default_output["runid"] = np.array([runid])
            default_output["ev"] = np.array([ev])
            # Check to see if there is a HumanGetBub file
            runidstr = str(runid[0])+"_"+str(runid[1])
            if runidstr != prev_runid:
                valid_file = False
                potential_bubble_file = os.path.join(bub_data_path, "HumanGetBub_"+str(runid[0])+"_"+str(runid[1])+".bin")
                if os.path.isfile(potential_bubble_file):
                    bub_data = RB(potential_bubble_file)
                    valid_file = True
                else:
                    to_add = deepcopy(default_output)
                    valid_file = False
            if valid_file:
                if np.count_nonzero(bub_data["ev"]==ev) == 0:
                    to_add = deepcopy(default_output)
                elif np.count_nonzero(bub_data["ev"]==ev) >= 1:
                    mask = np.where(bub_data["ev"]==ev)
                    to_add = {}
                    for k,v in bub_data.items():
                        to_add[k] = v[mask]
                else:
                    to_add = deepcopy(default_output)
            else:
                to_add = deepcopy(default_output)
            if first_file:
                out = to_add
                first_file = False
            else:
                out = dictionary_append(out, to_add)
            prev_runid = runidstr


        #
        #
        # for b_file in bubble_files:
        #     b_path = os.path.join(bub_data_path, "HumanGetBub_"+b_file+".bin")
        #     if first_file:
        #         out = RB(b_path)
        #         first_file = False
        #     else:
        #         out = dictionary_append(out, RB(b_path))
        # WB(output_path, out, rowdef=8)

    # debug = True
    # if debug:
    #     bubs = RB(output_path)
    #     timing = RB("/pnfs/coupp/persistent/grid_output/SBC-17/output/TimingAnalysis_all.bin")
    #
