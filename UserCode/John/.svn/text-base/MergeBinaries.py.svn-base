## John Gresl

import os
from collections import defaultdict, OrderedDict
import time
from copy import deepcopy

import numpy as np

from SBCcode.DataHandling.ReadBinary import ReadBlock as RB
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as WB

from SBCcode.AnalysisModules.AnalyzeDytran import dytranAnalysis as da
from SBCcode.AnalysisModules.EventAnalysis import EventAnalysis as eva
from SBCcode.AnalysisModules.ImageAnalysis import BubbleFinder
from SBCcode.AnalysisModules.AcousticT0 import AcousticAnalysis as aa
from SBCcode.AnalysisModules.PMTComprehensiveModule import PMTcm as pmtpa
from SBCcode.AnalysisModules.PMTfastDAQalignment import PMTandFastDAQalignment as pmtfda
from SBCcode.AnalysisModules.PTData import main as ptd
from SBCcode.AnalysisModules.TimingAnalysis import TimingAnalysis as ta
from SBCcode.DataHandling.GetSBCEvent import GetEvent as get_event
from SBCcode.DataHandling.WriteBinary import WriteBinaryNtupleFile as wb
from SBCcode.UserCode.John.NewT0 import calculate_t0 as calculate_t0


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
                print("DEBUG:", k, d1[k].shape, v.shape)
                d1[k] = np.append(d1[k], v, axis=0) # <-- If we have numpy arrays
    return d1


if __name__ == "__main__":
    file_templates = [#"AcousticAnalysis_{runid}.bin",
                      #"DytranAnalysis_{runid}.bin",
                      #"EventAnalysis_{runid}.bin",
                      #"HistoryAnalysis_{runid}.bin",
                      #"ImageAnalysis_{runid}.bin",
                      #"PMTfastDAQalignment_{runid}.bin",
                      #"PMTpulseAnalysis_{runid}.bin",
                      "TimingAnalysis_{runid}.bin",
                      ]
    defaults = [#aa(None, None),
                #da(None),
                #eva(None),
                #ptd(None),
                #BubbleFinder(None, None, None, None, None, None),
                #pmtfda(None),
                #pmtpa(None),
                ta(None, None, None),
                #calculate_t0(None, None, None, None)
                ]
    p_list = [(f, d) for f, d in zip(file_templates, defaults)]

    recon_directory =  "/pnfs/coupp/persistent/grid_output/SBC-17-T0Test3/output/"
    output_directory = "/nashome/j/jgresl/"
    runid_list = sort_runs([f for f in os.listdir(recon_directory) if os.path.isdir(os.path.join(recon_directory, f))])
    runid_list = trim_runlist(runid_list, start="20170623_3")
    bad_run_file = "BadRunsV6.npy"
    remake_badruns = False
    if remake_badruns:

        print("Building bad run list.")
        bad_runs = defaultdict()
        for f_temp in file_templates:
            bad_runs[f_temp] = set()
        bad_runs["AcousticAnalysis_{runid}.bin"] = set()
        for f_temp in file_templates:
            if f_temp != "AcousticTEST_{runid}.bin":
                for runid in runid_list:
                    if not os.path.isfile(os.path.join(recon_directory, runid, f_temp.format(runid=runid))):
                        print("\tSkipping {}. File not present."\
                              .format(os.path.join(recon_directory, runid, f_temp.format(runid=runid))))
                        bad_runs[f_temp].add(runid)
                        continue
                    try:
                        RB(os.path.join(recon_directory, runid, f_temp.format(runid=runid)), max_file_size=800)
                    except IndexError:
                        print("\tSkipping {}. File exists, but unable to read properly." \
                              .format(os.path.join(recon_directory, runid, f_temp.format(runid=runid))))
                        bad_runs[f_temp].add(runid)
                    except OSError:
                        print("\tSkipping {}. File above maximum file size. (Raise this in ReadBlock(...)" \
                              .format(os.path.join(recon_directory, runid, f_temp.format(runid=runid))))
                        bad_runs[f_temp].add(runid)
            if f_temp == "AcousticTEST_{runid}.bin":
                f_temp = "AcousticAnalysis_{runid}.bin"
                for runid in runid_list:
                    recon_directory = "/pnfs/coupp/persistent/grid_output/SBC-17-T0Test2/output"
                    if not os.path.isfile(os.path.join(recon_directory, runid, f_temp.format(runid=runid))):
                        print("\tSkipping {}. File not present."\
                              .format(os.path.join(recon_directory, runid, f_temp.format(runid=runid))))
                        bad_runs[f_temp].add(runid)
                        continue
                    try:
                        RB(os.path.join(recon_directory, runid, f_temp.format(runid=runid)), max_file_size=800)
                    except IndexError:
                        print("\tSkipping {}. File exists, but unable to read properly." \
                              .format(os.path.join(recon_directory, runid, f_temp.format(runid=runid))))
                        bad_runs[f_temp].add(runid)
                    except OSError:
                        print("\tSkipping {}. File above maximum file size. (Raise this in ReadBlock(...)" \
                              .format(os.path.join(recon_directory, runid, f_temp.format(runid=runid))))
                        bad_runs[f_temp].add(runid)
        print("----------")
        print("-Bad Runs-")
        print("----------")
        for k,v in bad_runs.items():
            print("{} has {} bad runs.".format(k, len(v)))
        print("\tSaving bad runs to {}.".format(bad_run_file))
        np.save(bad_run_file, bad_runs)
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ###############################################################################
    do_merge = True
    if do_merge:
        print("Beginning to merge all files.")
        big_out = []
        tstart = time.time()
        try:
            bad_runs = np.load(bad_run_file).flat[0]
        except:
            print("Couldn't find bad run file...")
            bad_runs = defaultdict(list)
        print("----------")
        print("-Bad Runs-")
        print("----------")
        for k, v in bad_runs.items():
            print("\t{} has {} bad runs.".format(k, len(v)))
        bad_list_intersection = set.intersection(*bad_runs.values())
        print("\t\tIntersection of bad events (failed for all analyses): {} items: {}".format(len(bad_list_intersection),
                                                                                              bad_list_intersection))
        for f_temp, d_default in p_list:
            if f_temp == "AcousticTEST_{runid}.bin":
                recon_directory = "/pnfs/coupp/persistent/grid_output/SBC-17-T0Test/output/"
                f_temp = "AcousticAnalysis_{runid}.bin"
            t0 = time.time()
            out = {}
            first_real_run = True
            print("\tStarting {}".format(f_temp))
            n_to_process = len(runid_list) - len(bad_runs[f_temp])
            n_processed = 0
            for n, runid in enumerate(runid_list):
                npev = np.array([int(runid.split("_")[1])], dtype=np.int32)
                nprunid = np.int32(runid.split("_"))
                print(nprunid)
                if n%25 == 0: # Print a status message every 25 runs
                    print("\t\t{:.2f}% done with file {}".\
                          format(n_processed/n_to_process*100, f_temp))
                if runid in bad_list_intersection: # If it failed for all analysis, skip it entirely.
                    continue
                if runid in bad_runs[f_temp]:
                    # print("#####FAKE VALUES#####")
                    temp = deepcopy(d_default)
                    temp["runid"] = np.array([nprunid])
                    temp["ev"] = npev
                    if first_real_run:
                        out = temp
                        first_real_run = False
                        # big_out.append(out)
                    else:
                        # big_out.append(temp)
                        out = dictionary_append(out, temp)
                    continue
                # print("\t\t\tRunID: {}".format(runid))

                if first_real_run:
                    out = RB(os.path.join(recon_directory, runid, f_temp.format(runid=runid)), max_file_size=800)
                    # big_out.append(out)
                    first_real_run = False
                else:
                    # print("******REAL VALUES*******")
                    d=RB(os.path.join(recon_directory, runid, f_temp.format(runid=runid)),
                                                    max_file_size=800)
                    # big_out.append(d)
                    try:
                        out = dictionary_append(out, d)
                    except:
                        print("Failed.#########")
                        if runid == "20170805_4":
                            for eee in range(12):
                                temp = deepcopy(d_default)
                                temp["runid"] = np.array([nprunid])
                                temp["ev"] = npev
                                out = dictionary_append(out, temp)
                        else:
                            temp = deepcopy(d_default)
                            temp["runid"] = np.array([nprunid])
                            temp["ev"] = npev
                            out = dictionary_append(out, temp)
                n_processed += 1

            t1 = time.time()
            print("\tTook {:.2f} seconds to read input files for {}".format(t1 - t0, f_temp))
            print("\tStarting to write output merged file.")
            WB(os.path.join(output_directory, f_temp.format(runid="all")), out, rowdef=1)
            t2 = time.time()
            print("\tTook {:.2f} seconds to write merged file {}".format(t2 - t1, f_temp.format(runid="all")))
            print("\tTook {:.2f} seconds for entire process to read and create {}". \
                  format(t2 - t0, f_temp.format(runid="all")))
        tfinish = time.time()
