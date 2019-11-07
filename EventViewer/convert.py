from glob import glob
import numpy as np
import os
import re
import time

skip = ['timestamp', 'livetime', 'piezo_max(3)', 'piezo_min(3)', 'piezo_starttime(3)', 'piezo_endtime(3)', 'piezo_freq_binedges(9)', 'acoustic_neutron', 'acoustic_alpha', 'scanner_array(2)', 'scan_source_array(2)', 'scan_nbub_array(2)', 'scan_trigger_array(2)', 'scan_comment_array(2)', 'scaler(8)', 'led_max_amp(8)', 'led_max_time(8)', 'null_max_amp(8)', 'first_hit(8)', 'last_hit(8)', 'max_amps(8)', 'max_times(8)', 'nearest_amps(8)', 'nearest_times(8)', 'numtrigs(8)', 'numpretrigs(8)', 'scan_comment_array(2)']
dtypes = {'s': 'U12', 'd': 'i4', 'f': 'f4', 'e': 'f4'}  # map fscanf format to numpy datatype
reco_directory = './'
raw_directory = './'
reco_directory = '/bluearc/storage/recon/current/30l-16/output'
raw_directory = '/bluearc/storage/30l-16-data/'

def natural_sort(things):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(things, key=alphanum_key)

def load_reco(filename):
    path = os.path.join(reco_directory, filename)
    with open(path) as f:
        files = f.readline()
        fields = f.readline().split()
        formats = f.readline().split()

        dt = []
        columns = []
        column = 0
        for field in fields:
            format = formats[column]
            for x in dtypes:
                if x in format:
                    dtype = dtypes[x]
            match = re.search('.*\((.*)\)', field)
            if match:
                dimensions = [int(x) for x in match.groups()[0].split(',')]
                length = np.prod(dimensions)
                types = formats[column:column + length]
                if len(set(types)) > 1:
                    message = 'cannot parse {} with types {} because mixed types are not supported'
                    raise NotImplementedError(message.format(field, types))
                if (len(dimensions) == 1) and (field not in skip):  # skip loading multidimensional arrays to save memory
                    columns.extend(list(range(column, column + length)))
                    dt.append((field, (dtype, dimensions)))
                column += length
            else:
                if field not in skip:
                    dt.append((field, dtype))
                    columns.append(column)
                column += 1
    try:
        old_events = np.load(os.path.join(reco_directory, 'reco_events.npy'))
    except FileNotFoundError:
        old_events = np.empty((0,), dtype=np.dtype(dt))

    skip_header = 6 + len(old_events)
    print('skipping {} reco lines which have already been parsed'.format(len(old_events)))
    new_events = np.genfromtxt(path, dtype=old_events.dtype, delimiter="  ", skip_header=skip_header, usecols=columns)

    return np.concatenate([old_events, new_events])

def validate(events):
    res = []
    for run, event, index in events:
        path = os.path.join(raw_directory, run, str(event), 'Event.txt')
        if os.path.isfile(path):
            res.append((run, event, index))
        else:
            print('not found {}'.format(path))

    return np.array(res, dtype=events.dtype)

def load_raw(filename, reco):
    try:
        old_events = np.load(os.path.join(reco_directory, filename))
    except FileNotFoundError:
        old_events = np.array([], dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])

    events = []
    runs = [os.path.basename(x) for x in glob(os.path.join(raw_directory, "20*"))]
    for run in runs:
        for event in natural_sort(glob(os.path.join(raw_directory, run, '[0-9]*/'))):
            event = os.path.basename(event.strip('/'))
            matches = np.argwhere((reco['run'] == run) & (reco['ev'] == int(event))).flatten()
            index = matches[0] if len(matches) > 0 else -1
            events.append((run, event, index))

    new_events = np.setdiff1d(np.array(events, dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')]), old_events)
    print('new events to be added: {}'.format(len(new_events)))
    new_events = validate(new_events)

    return np.concatenate([old_events, new_events])


print('starting now')
start = time.time()
reco = load_reco('merged_all.txt')
np.save(os.path.join(reco_directory, 'reco_events'), reco)

raw = load_raw('raw_events.npy', reco)
np.save(os.path.join(reco_directory, 'raw_events'), raw)
print('finished in {:.0f} seconds'.format(time.time() - start))
