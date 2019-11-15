# John Gresl

import sys
import os
import time
import re
import numpy as np


import SBCcode.DataHandling.ReadBinary

skip = ["timestamp", "livetime", "piezo_max(3)", "piezo_min(3)",
        "piezo_starttime(3)", "piezo_endtime(3)",
        "piezo_freq_binedges(9)", "acoustic_neutron",
        "acoustic_alpha", "scanner_array(2)", "scan_source_array(2)",
        "scan_nbub_array(2)", "scan_trigger_array(2)",
        "scan_comment_array(2)", "scaler(8)", "led_max_amp(8)",
        "led_max_time(8)", "null_max_amp(8)", "first_hit(8)",
        "last_hit(8)", "max_amps(8)", "max_times(8)",
        "nearest_amps(8)", "nearest_times(8)", "numtrigs(8)",
        "numpretrigs(8)", "piezo_time_windows", "piezoE",
        "PressureBins", "TempData", "PMTmatch_area_nobs",
        "nVetohits_fastdaq", "nPMThits_fastdaq", "PMTmatch_min",
        "PMTmatch_max", "PMTmatch_area"]

dtypes = {'s': 'U12', 'd': 'i4', 'f': 'f4', 'e': 'f4',
          "uint32": "u4", "int32": "i4", "float64": "f8",
          "int64": "i8", "int16": "i2", "int8": "i", "<U32": "U12"}  # map fscanf format to numpy datatype


def _compare(names, types):
    out_str = ""
    format_str = "{{:{}s}}"
    lens = []
    for name, typ in zip(names, types):
        lens.append(str(max([len(str(name)), len(str(typ))])))
    out_str += "Names: | "
    for name, lenx in zip(names, lens):
        out_str += format_str.format(lenx).format(name)+" | "
    out_str += "\nTypes: | "
    for typ, lenx in zip(types, lens):
        out_str += format_str.format(lenx).format(str(typ))+" | "
    out_str += "\n"
    return out_str


def np_to_py(var):
    if type(var) == np.ndarray:
        return var.tolist()
    if type(var) == np.float64:
        return float(var)
    if type(var) == np.int32:
        return int(var)
    if type(var) == np.int64:
        return int(var)
    if type(var) == np.str_:
        return str(var)
    raise TypeError("Type {} not implemented.".format(type(var)))


def bin_to_npy(reco_bin_file, np_save_output, verbose=False):
    def _verbose(s):
        if verbose:
            print(s)
        return
    _verbose("Using file '{}' to create '{}'".format(reco_bin_file,
                                                     np_save_output))
    start_time = time.time()
    data = SBCcode.DataHandling.ReadBinary.ReadBlock(reco_bin_file)
    read_time = time.time()
    _verbose("Time to read the file {}: {}s".format(input, read_time-start_time))
    temp_keys = list(data.keys())
    for k in temp_keys:
        if k in skip:
            _verbose("Deleting key {}".format(k))
            del data[k]
    del temp_keys
    data_keys = list(data.keys())

    n_events = len(data[data_keys[0]])
    temp_runid2 = np.zeros(n_events).astype("U")
    if "runid" in data_keys:
        for i in range(n_events):
            temp_runid2[i] = "_".join([str(data["runid"][i][j]) for j in range(2)])
        data["runid"] = temp_runid2
    value_lengths = []
    new_dtypes = []
    for key in data_keys:
        value_lengths.append(data[key].shape[0])
        if len(data[key].shape) == 1:
            _verbose("Key {} with value {} is being cast as type {}".
                     format(key, data[key][0], dtypes[str(data[key].dtype)]))
            new_dtypes.append((key, dtypes[str(data[key].dtype)]))
        elif len(data[key].shape) == 2:
            _verbose("Key {} with value {} is being cast as type {}".
                     format(key, data[key][0], (dtypes[str(data[key].dtype)], [data[key].shape[1]])))
            new_dtypes.append((key, (dtypes[str(data[key].dtype)], [data[key].shape[1]])))
        else:
            _verbose(str(key)+" needs to be looked at...")
            print(data[key].shape)
            raise ValueError()
    if not all([value_lengths[0] == val for val in value_lengths]):
        raise ValueError("# of events is not consitent across each key." +
                         "No defined behaviour for this. *ABORTING*")
    n_events = value_lengths[0]
    out = []
    dts = []
    for n in range(n_events):
        tmp = []
        dts = []
        for (k, v), t in zip(list(data.items()), new_dtypes):
            tmp.append(np_to_py(v[n]))
            dts.append(t)
        tmp = tuple(tmp)
        out.append(tmp)
    _verbose(_compare(list(data_keys), new_dtypes))
    out = np.array(out, dtype=np.dtype(dts))
    np.save(np_save_output, out)
    process_time = time.time()
    _verbose("Time to process {} events: {}s".format(n_events, process_time-start_time))
    return


def txt_to_npy(reco_txt_file, np_save_output, verbose=False):
    def _verbose(s):
        if verbose:
            print(s)
        return
    _verbose("Using file '{}' to create '{}'".format(reco_txt_file,
                                                     np_save_output))
    start_time = time.time()
    with open(reco_txt_file) as f:
        files = f.readline()  # <-- Unused but it pushes the buffer forward.. and it IS the string of filenames anyway.
        fields = f.readline().split()
        formats = f.readline().split()
        dt = []
        columns = []
        cur_column = 0
        for field in fields:
            form = formats[cur_column]
            for x in dtypes:
                if x in form:
                    dtype = dtypes[x]
            match = re.search('.*\((.*)\)', field)
            if match:
                dimensions = [int(x) for x in match.groups()[0].split(',')]
                length = np.prod(dimensions)
                types = formats[cur_column:cur_column + length]  # <-- Not important
                if len(set(types)) > 1:
                    print("Cannot parse {} with types {} because mixed types are not supported".format(field, types))
                    return
                if (len(dimensions) == 1) and (field not in skip):  # skip loading nd-arrays to save memory
                    columns.extend(list(range(cur_column, cur_column + length)))
                    dt.append((field, (dtype, dimensions)))
                cur_column += length
            else:
                if field not in skip:
                    dt.append((field, dtype))
                    columns.append(cur_column)
                cur_column += 1
    old_events = np.empty((0,), dtype=np.dtype(dt))
    skip_header = 6 + len(old_events)
    _verbose('skipping {} reco lines which have already been parsed'.format(len(old_events)))
    out = np.genfromtxt(reco_txt_file, dtype=old_events.dtype, delimiter="  ", skip_header=skip_header,
                        usecols=columns)
    np.save(np_save_output, out)
    process_time = time.time()
    _verbose("Successfully processed {} events in {}s".format(len(out), process_time-start_time))
    return


def main(*args, mode="binary", use_shell=True, verbose=False):
    if use_shell:
        if len(args) != 3:
            print("Incorrect usage. Proper usage is:\n")
            print("python bin_to_npy.py <directory-of-merged_all-file> -mode=<text or binary>")
            return
        reco_dir = args[1]
        print(args)
        mode = args[2].split("=")[-1]
        if mode not in ["text", "binary"]:
            print("Invalid mode: '{}'. Must be either 'text' or 'binary'".format(mode))
            return
    else:
        reco_dir = args[0]
    file_in = os.path.join(reco_dir, "merged_all." + "bin" if mode == "binary" else "txt")
    npy_out = os.path.join(reco_dir, "reco_events.npy")
    if not os.path.isfile(file_in):
        print("No file exists at '{}'\nAborting.".format(file_in))
        return
    if os.path.isfile(npy_out):
        while True:
            proceed = input("\nWARNING. A file already exists at '{}'\n".format(npy_out) +
                            "This procedure will overwrite the file.\n" +
                            "Proceed? (y/n): ").lower()
            if proceed == "y":
                break
            if proceed == "n":
                print("Aborting.")
                return
            print("Invalid input: '{}'".format(proceed))
    if not os.path.isdir(os.path.dirname(npy_out)):
        print("Invalid path of output file: '{}'".format(npy_out))
        return
    if mode == "binary":
        bin_to_npy(file_in, npy_out, verbose=verbose)
    if mode == "text":
        txt_to_npy(file_in, npy_out, verbose=verbose)
    return

if __name__ == "__main__":
    # If you would like to hard-code the paths in this file, you can! Set the use_command_line
    # variable to False and then set the bin_in and npy_out paths appropriately.
    use_command_line = True     # If use_command_line is True, then bin_in and npy_out will
                                # automatically be replaced with the command line arguments.
    verbose = True
    mode = "binary"
    reco_dir = "/pnfs/coupp/persistent/grid_output/SBC-17/output"
    npy_out = "/pnfs/coupp/persistent/grid_output/SBC-17/output/reco_events.npy"
    if use_command_line:
        main(*sys.argv, use_shell=True, verbose=verbose)
    else:
        main(reco_dir, mode=mode, use_shell=False, verbose=verbose)
