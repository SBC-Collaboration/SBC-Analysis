# This version only creates the raw npy file, but takes as input the reco if available, in order to better sort the events

# Usage example: python convert.py /bluearc/storage/recon/current/30l-16/output /bluearc/storage/30l-16-data
# Produces npy files for navigation from raw data, to be used by PED event display
# may need to source /coupp/data/home/coupp/PEDsvn/setup_ped_paths.sh

from glob import glob
import numpy as np
import os
import re
import time
import sys
import logging

skip = ['timestamp', 'livetime', 'piezo_max(3)', 'piezo_min(3)', 'piezo_starttime(3)', 'piezo_endtime(3)', 'piezo_freq_binedges(9)', 'acoustic_neutron', 'acoustic_alpha', 'scanner_array(2)', 'scan_source_array(2)', 'scan_nbub_array(2)', 'scan_trigger_array(2)', 'scan_comment_array(2)', 'scaler(8)', 'led_max_amp(8)', 'led_max_time(8)', 'null_max_amp(8)', 'first_hit(8)', 'last_hit(8)', 'max_amps(8)', 'max_times(8)', 'nearest_amps(8)', 'nearest_times(8)', 'numtrigs(8)', 'numpretrigs(8)', 'scan_comment_array(2)']
dtypes = {'s': 'U12', 'd': 'i4', 'f': 'f4', 'e': 'f4'}  # map fscanf format to numpy datatype

if len(sys.argv) != 2:
    print('Should be 1 arguments.')
    print('Usage: python convert.py <dir-where-raw-data-is>')
    exit()

# print('Number of arguments:' + str(len(sys.argv)) + 'arguments.')
# print('Argument List:' + str(sys.argv))

raw_directory = str(sys.argv[1])

# print('reco_directory = ' + reco_directory)
# print('raw_directory = ' + raw_directory)

def natural_sort(things):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(things, key=alphanum_key)

# In this version, don't create a new reco, just load it, to use in load_raw if available
def load_reco(filename):
    try:
        old_events = np.load(os.path.join(reco_directory, 'reco_events.npy'))
    except FileNotFoundError:
        old_events = np.empty((0,), dtype=np.dtype(dt))

    return old_events

def validate(events):
    res = []
    for run, event, index in events:
        path = os.path.join(raw_directory, run, str(event), 'Event.txt')
        if os.path.isfile(path):
            res.append((run, event, index))
        else:
            print('not found {}'.format(path))

    return np.array(res, dtype=events.dtype)

def load_raw(filename, reco=None):
    try:
        old_events = np.load(os.path.join(raw_directory, filename))
    except FileNotFoundError:
        old_events = np.array([], dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')])

    events = []
    runs = [os.path.basename(x) for x in glob(os.path.join(raw_directory, "20*"))]
    for run in runs:
        for event in natural_sort(glob(os.path.join(raw_directory, run, '[0-9]*/'))):
            event = os.path.basename(event.strip(os.sep))
            if reco is not None:
                matches = np.argwhere((reco['run'] == run) & (reco['ev'] == int(event))).flatten()
                index = matches[0] if len(matches) > 0 else -1
            else:
                index = -1
            events.append((run, event, index))

    new_events = np.setdiff1d(np.array(events, dtype=[('run', 'U12'), ('ev', 'i4'), ('reco index', 'i4')]), old_events)
    print('new events to be added: {}'.format(len(new_events)))
    new_events = validate(new_events)

    return np.concatenate([old_events, new_events])


print('starting now')
start = time.time()

try:
    raw = load_raw('raw_events.npy')
    np.save(os.path.join('/nashome/z/zsheng/Event-Viewer/raw_events.npy'), raw)
except:
    print('Cannot find raw data, skipping creation of raw_events.npy\n')

print('finished in {:.0f} seconds'.format(time.time() - start))
