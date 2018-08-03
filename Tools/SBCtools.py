# This file is meant to be used as a set of utilities useful for others. Please try to follow
# the example format below. It's probably a good idea to raise exceptions in this file, and then
# decide how you want to handle them in your own scripts. If possible, provide tests immediately after
# in an if __name__ == "__main__": block. Include normal and edge cases.


def example_function(s, n, sep=""):
    # Author: John Gresl 8/2/2018
    # Inputs:
    #   s: Any string.
    #   n: An integer. Number of times to copy s.
    #   sep: A string. Will be inserted between every copy of s.
    # Outputs:
    #   The string s repeated n times.
    # Example: example_function("hello", 3) == "hellohellohello"
    return ((s + sep) * n)[:-len(sep)]


if __name__ == "__main__":
    assert example_function("hello", 3, sep=" ") == "hello hello hello"
    assert example_function("", 5, sep="x") == "xxxx"
    assert example_function("hello", 0) == ""
    assert example_function("dark matter", 2, sep=" woah ") == "dark matter woah dark matter"

"""
Table of Contents:

real_number_q(n): Returns true if n is a number (numpy int, float, uint, or python int or float)

dictionary_append(d1, *args, make_copy=True): Returns a dictionary with the values from dictionaries in args
                                              appended to the values in d1. Modifies d1 in place if make_copy is False,
                                              otherwise returns a deepcopy of d1.

sort_runs(arr, reverse): Sorts an array of runs if the array looks like ["20170623_0", "20170623_5", ...]

get_runs(dir, search_for="folders", do_sort=True, reverse=False): Compiles a list of runids from a directory.

trim_run_list(arr, start=None, stop=None): Trims a runlist like one that was returned from get_runs or sort_runs.
                                           start/stop should be a runid like "20170623_0" or None.


"""

# Native python imports
import os
from copy import deepcopy
from collections import defaultdict, OrderedDict
import operator
from random import shuffle
import re

# Anaconda python imports
import numpy as np


# Custom python imports


def real_number_q(n):
    # Author: John Gresl 8/2/2018
    # Inputs:
    #   n: Any variable
    # Outputs: A boolean value. True if n is any type of number (numpy or python versions). False otherwise
    return type(n) in [int, float] or issubclass(type(n), (np.int8, np.int16, np.int32, np.int64,
                                                           np.uint8, np.uint16, np.uint32, np.uint64,
                                                           np.float16, np.float32, np.float64))
if __name__ == "__main__":
    assert real_number_q(0)
    assert real_number_q(1.2)
    assert real_number_q(-5)
    assert real_number_q(np.int32(5))
    assert real_number_q(np.uint32(5123))
    assert real_number_q(np.float16(4))
    assert real_number_q(np.uint64(2147483647))
    assert not real_number_q("5")
    assert not real_number_q("dark matter")


def dictionary_append(d1, *args, make_copy=True):
    # Author: John Gresl 8/2/2018
    # Inputs:
    #   d1: Dictionary
    #   args: Any number of dictionaries
    #   make_copy: If this flag is set to True, this will create a new dictionary. If this flag is false,
    #              then d1 will be modified
    # Outputs: A dictionary with the same keys, but the values are the values from args appended to the values of d1
    # Note: The keys in d1 and args MUST match, and all of the values MUST be lists.
    if make_copy:
        d1 = deepcopy(d1)
    if len(args) == 0:
        raise TypeError("dictionary_append must be called with at least 2 arguments.")
    for arg in args:
        if type(arg) not in [dict, defaultdict, OrderedDict]:
            raise TypeError("args must be dictionaries!")
    for d2 in args:
        if not set(d1.keys()) == set(d2.keys()):
            raise KeyError("The keys for the two dictionaries must match! Mismatched keys = {}". \
                           format(set(d1.keys()).symmetric_difference(set(d2.keys()))))
    for d2 in args:  # This '2nd' for loop because we want to make sure all the keys are the same first.
        for k, v in d2.items():
            if type(v) not in [list, np.ndarray]:
                raise ValueError("The values of the dictionary MUST be list or np.ndarray. Key {} has type(value)={}". \
                                 format(k, type(v)))
            try:
                d1[k].extend(v)  # <-- If we have python lists as values
            except AttributeError:
                # print("DEBUG:", k, d1[k].shape, v.shape) # Uncomment this line to help you debug your dictionaries
                d1[k] = np.append(d1[k], v, axis=0)  # <-- If we have numpy arrays
    return d1
if __name__ == "__main__":
    d1 = {"a": [1, 2, 3], "b": [6, 1, 2], "c": [0, 0, 1]}
    d2 = {"a": [6, 2], "b": [1, 2], "c": [9, 2]}
    d3 = {"a": [], "b": [], "c": []}
    d4 = {"a": [1, 2], "b": [6, 2], "c": [6, 1], "x": [9, 9]}
    assert dictionary_append(d1, d2, make_copy=True) == {"a": [1, 2, 3, 6, 2], "b": [6, 1, 2, 1, 2],
                                                         "c": [0, 0, 1, 9, 2]}
    assert d1 is not dictionary_append(d1, d2, make_copy=True) == {"a": [1, 2, 3, 6, 2], "b": [6, 1, 2, 1, 2],
                                                                   "c": [0, 0, 1, 9, 2]}
    assert dictionary_append(d1, d3, make_copy=True) == d1
    assert d1 is dictionary_append(d1, d3, make_copy=False)
    try:
        dictionary_append(d1, d4, make_copy=True)
    except KeyError:
        pass


def sort_runs(arr, reverse=False):
    # Author: John Gresl 8/2/2018
    # Input:
    #   arr: An array of run_ids as strings. Should look like ["20170623_0", "20170623_5", etc...]
    #   reverse: Bool. If true, returns the reverse-sorted list.
    # Outputs: A natural-ish sorted version that puts the dates in order and the run numbers for each date in order
    s = sorted([np.int32(runid.split("_")) for runid in arr], key=operator.itemgetter(0, 1), reverse=reverse)
    return ["_".join(np.str(runid).strip("[]").split()) for runid in s]
if __name__ == "__main__":
    test_array = []
    for i in range(20):
        for j in range(20):
            test_array.append("201606" + str(10 + i) + "_" + str(j))
            # The list is created already sorted.
    assert sort_runs(test_array) == test_array
    test_array_copy = deepcopy(test_array)
    shuffle(test_array_copy)  # random.shuffle modifies the list in place.
    assert sort_runs(test_array_copy) == test_array


def get_runs(dir, search_for="folders", do_sort=True, reverse=False):
    # Author: John Gresl
    # Inputs:
    #   dir: A directory containing folders of all runs
    #   search_for: A string indicating whether to extract run_ids from folders or from files
    #               Can be either "folders" or "files"
    #   do_sort: A boolean. If true, will return the sorted run_list from sort_runs above.
    #   reverse: Bool. If true, returns the reverse-sorted list.
    # Outputs: An array of run-strings from a given directory. This is passed to sort_runs if do_sort is true
    if search_for.lower().strip() == "folders":
        base_dirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    elif search_for.lower().strip() == "files":
        base_dirs = [d for d in os.listdir(dir) if os.path.isfile(os.path.join(dir, d))]
    else:
        raise ValueError("'search_for' must be either 'folders' or 'files. Got {}".format(search_for))
    to_match = re.compile("([0-9]{8}_[0-9]+){1}")
    out = []
    for d in base_dirs:
        match_q = re.search(to_match, d)
        if match_q:
            out.append(match_q.group(0))
    if do_sort:
        return sort_runs(out, reverse=reverse)
    return out
if __name__ == "__main__":
    file_dir = "/coupp/data/home/coupp/HumanGetBub_output_SBC-17/"
    folder_dir = "/pnfs/coupp/persistent/grid_output/SBC-17/output"
    print(get_runs(file_dir, search_for="files", do_sort=True))
    print(get_runs(folder_dir, search_for="folders", do_sort=True))



def trim_runlist(arr, start=None, stop=None):
    # Author: John Gresl 8/2/2018
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