## John Gresl



import numpy as np
from copy import deepcopy
import os
from collections import defaultdict
from collections import OrderedDict
import re

from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as WB
from SBCcode.DataHandling.GetSBCEvent import GetEvent as GE


def sort_runs(arr):
    # Input:
    #   arr: An array of run_ids as strings. Should look like ["20170623_0", "20170623_5", etc...]
    # Outputs: A natural-ish sorted version that puts the dates in order and the run numbers for each date in order
    dates_only = []
    runs_only = []
    for run_id in arr:
        dates_only.append(run_id.split("_")[0])
        runs_only.append(run_id.split("_")[1])
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


def trim_runlist(arr, start=None, stop=None):
    # Inputs:
    #   arr: An array of run_ids as strings. Should look like ["20170623_0", "20170623_5", etc...]
    #   start: Start run number. If this is not supplied, will start at the beginning
    #   stop: Stop run number. If this is not supplied, will continue to end
    # Outputs: A sorted, trimmed runlist that goes from start to stop
    arr = sort_runs(arr)
    start = arr[0] if start == None else start
    stop = arr[-1] if stop == None else stop
    start_date = int(start.split("_")[0])
    start_run_num = int(start.split("_")[1])
    stop_date = int(stop.split("_")[0])
    stop_run_num = int(stop.split("_")[1])

    out = [ ]
    for run in arr:
        date = int(run.split("_")[0])
        run_num = int(run.split("_")[1])
        if start_date > date or date > stop_date:
            continue
        if (start_run_num > run_num and date == start_date) or (run_num > stop_run_num and date == stop_date):
            continue
        out.append(run)
    return out

def dictionary_append(d1, *args):
    # Inputs:
    #   d1: Dictionary
    #   args: Any number of dictionaries
    # Outputs: A dictionary with the same keys, but the values are the values from args appended to the values of d1
    # Note: The keys in d1 and args MUST match, and all of the values MUST be lists.
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
                d1[k] = np.append(d1[k], v, axis=0) # <-- If we have numpy arrays
    return d1


if __name__ == "__main__":
    ### To anyone (and myself) who uses this script in the future!!!!!!::::
    ### I noticed that there are some events in the HumanGetBub files that did NOT appear in the timing analysis
    ### files. (There might be a better way to do this, but this is what worked for me):
    ### To make sure all the dimensions line up when using the cut, you can use the commented out code at the bottom
    ### and look at the final output where it says "The difference in elements is: " and add those elements to the
    ### BAD_DICT variable. The key is the runid, and the value is a list of bad events for that runid.

    RAW_DIRECTORY = "/bluearc/storage/SBC-17-data/"
    BUBBLE_DIRECTORY = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/"
    RECON_DIRECTORY = "/pnfs/coupp/persistent/grid_output/SBC-17/output/"
    ma = RB("/pnfs/coupp/persistent/grid_output/SBC-17/output/TimingAnalysis_all.bin")
    START_DATE = "20170919"
    START_RUN = START_DATE+"_0"
    END_DATE = None
    BAD_LIST = ['20170923_2_45', '20170929_0_24', '20170925_4_0', '20171009_3_0',
                '20170929_3_37', '20171007_3_11', '20171008_5_42', '20170925_3_0',
                '20171008_2_26', '20170924_3_79', '20170925_2_50', '20170930_6_0',
                '20170927_1_43', '20171008_3_0', '20170920_5_2']


    BAD_DICT = defaultdict(list,
                           **{"20170923_2": [45], #
                            "20170929_0": [24], #
                #"20170925_4": [0],
                #"20171009_3": [0],
                "20170929_3": [37], #
                "20171007_3": [11], #
                "20171008_5": [42], # asdasd
                #"20170925_3": [0],
                "20171008_2": [26], #
                "20170924_3": [79], #
                "20170925_2": [50], #
                #"20170930_6": [0],
                "20170927_1": [43], #
                #"20171008_3": [0],
                "20170930_6": [1],  #
                "20170920_5": [2]}) #
    #BAD_DICT = defaultdict(list)
    BAD_EVS = []
    #BAD_LIST = []

    bub_files = [ f for f in os.listdir(BUBBLE_DIRECTORY) if os.path.isfile(f)]
    runids = []
    for bub_file in bub_files:
        runids.append(bub_file.split("_")[1]+"_"+bub_file.split("_")[2])


    # 1. Build a sorted list of events.
    sub_drs = [d for d in os.listdir(RECON_DIRECTORY) \
               if os.path.isdir(os.path.join(RECON_DIRECTORY, d))]
    sorted_drs = sort_runs(sub_drs)
    my_run_list = trim_runlist(sorted_drs, start=START_RUN, stop=None)
    my_runlist = []
    for r1 in ma["runid"]:
        my_runlist.append(str(r1[0])+"_"+str(r1[1]))
    my_run_list = sort_runs(list(set(my_runlist)))
    print(my_run_list)
    first_event = True
    for event in my_run_list:
        #event_path = os.path.join(RECON_DIRECTORY, event)
        event_path = BUBBLE_DIRECTORY
        data_path = os.path.join(event_path, "HumanGetBub_"+event+".bin")
        raw_data_path = os.path.join(RAW_DIRECTORY, event)
        # We need to know how many events there were for the default output, so we look in the raw directory
        # As of writing this, an event is stored denoted by a folder with an integer name. Change the
        # following match if this changes in the future. Must also change the default output in case you change
        # to using non-ints as event number specifiers
        event_re = re.compile("^\d+$") # <-- Should match any string of 1 or more digits
        matched_events = [mtch for mtch in os.listdir(raw_data_path) if \
                          (event_re.match(mtch) and (int(mtch) not in BAD_DICT[event]))] # if event in BAD_DICT.keys() else True))]
        matched_events = []
        skip = False
        for mtch in os.listdir(raw_data_path):

            if event_re.match(mtch):
                if int(mtch) not in BAD_DICT[event]:
                    matched_events.append(mtch)
                else:
                    print("\nskipping event {} in run {}".format(mtch, event))
                    skip = True

        if skip:
                print(matched_events)
        # print("Matched events for {}: {}".format(event, matched_events))
        # print("Bad Events: {}".format(BAD_DICT[event]))
        # print()
        n_events = len(matched_events)

        default_output = OrderedDict()
        default_output["runid"]=\
            np.repeat([np.array([int(event.split("_")[0]), int(event.split("_")[1])], dtype=np.uint32)], n_events, axis=0)
        default_output["ev"]= np.array(sorted(matched_events, key=int), dtype=np.int32)
        if skip:
            print(default_output["ev"])
            print()
        default_output["ibubimage"]= np.repeat(np.int32(-1), n_events)
        default_output["nbubimage"]= np.repeat(np.int32(-1), n_events)
        default_output["cam"]= np.repeat(np.int32(-1), n_events)
        default_output["frame"]= np.repeat(np.int32(-1), n_events)
        default_output["ipix"]= np.repeat(np.int32(-1), n_events)
        default_output["jpix"]= np.repeat(np.int32(-1), n_events)
        run_out = deepcopy(default_output)
        total_out = default_output
        try:
            bub_data = RB(data_path)
            if skip:
                matched_events = [int(i) for i in matched_events]
                dif = set(matched_events).symmetric_difference(set(bub_data["ev"]))
                for bad_event in dif:
                    for n in range(len(bub_data["runid"])):
                        if str(bub_data["runid"][n][0])+"_"+str(bub_data["runid"][n][1]) == event and bub_data["ev"][n] == bad_event:
                            for k in bub_data.keys():
                                bub_data[k] = np.delete(bub_data[k],n, axis=0)

            if first_event:
                total_out = run_out
                first_event = False
            else:
                total_out = dictionary_append(total_out, bub_data)
                #print("DEBUG2:", total_out["runid"])

        except FileNotFoundError as e:
            print("Unable to load {}. Appending default value.".format(event))
            total_out = dictionary_append(total_out, default_output)
        except Exception as e:
            print("Uh oh? Something really bad went wrong?")
            raise e
    WB(os.path.join(RECON_DIRECTORY, "HumanGetBub_all.bin"), total_out)
    ################################################################################################
    ################################################################################################
    ################################################################################################
    ################################################################################################
    ################################################################################################
    ################################################################################################
    ################################################################################################
    ################################################################################################
    ## Debugging
    ma = RB("/pnfs/coupp/persistent/grid_output/SBC-17/output/TimingAnalysis_all.bin")
    ma_arr  = []
    for r,e in zip(ma["runid"], ma["ev"]):
        ma_arr.append(str(r[0])+"_"+str(r[1])+"_"+str(e))

    to_arr = []
    for r,e in zip(total_out["runid"], total_out["ev"]):
        to_arr.append(str(r[0])+"_"+str(r[1])+"_"+str(e))

    ma_set = set(ma_arr)
    to_set = set(to_arr)

    print("Length of ta_arr is {}".format(len(ma_arr)))
    print("Length of to_arr is {}".format(len(to_arr)))
    ma_dict = OrderedDict()
    to_dict = OrderedDict()
    help_me = 0
    for e in ma_set:
        ma_dict[e] = ma_arr.count(e)
        to_dict[e] = to_arr.count(e)
        if ma_arr.count(e) != to_arr.count(e):
            #print("Merged all -- {} -- {}".format(e, ma_arr.count(e)))
            #print("Total out  -- {} -- {}".format(e, to_arr.count(e)))
            #print()
            help_me += abs(ma_arr.count(e) - to_arr.count(e))
    print("Difference in length: {}".format(len(ma_arr) - len(to_arr)))
    print("Other difference in length: {}".format(help_me))

    print("Length of ta_set is {}".format(len(ma_set)))
    print("Length of to_set is {}".format(len(to_set)))
    print("The differences in elements is {}".format(ma_set.symmetric_difference(to_set)))

