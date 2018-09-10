import os
from .ReadBinary import ReadBlock as rb
from .ReadText import ReadFile as rt
import traceback
import numpy as np


full_loadlist = ['fastDAQ', 'slowDAQ', 'PMTtraces', 'event',
                 'camdata', 'images', 'DAQsetings']


def GetEvent(rundirectory, ev, *loadlist, max_file_size=None):
    event = dict()
    event_dir = os.path.join(rundirectory, str(ev))
    for key in full_loadlist:
        event[key] = dict(loaded=False)

    neglist = False
    if len(loadlist) == 0:
        loadlist = full_loadlist
    elif len(loadlist[0]) > 0 and loadlist[0][0:1] == '~':
        neglist = True

    if ('fastDAQ' in loadlist) or (neglist and '~fastDAQ' not in loadlist):
        i_file = 0
        while True:

            binfile = os.path.join(event_dir,
                                   'fastDAQ_' + str(i_file) + '.bin')
            calfile = os.path.join(event_dir,
                                   'fastDAQ_' + str(i_file) + '_cal.txt')

            if not (os.path.exists(binfile) and os.path.exists(calfile)):
                break

            try:
                d_bin = rb(binfile, max_file_size=max_file_size)
                d_cal = rt(calfile)
                d = dict()
                for key in d_bin:
                    d[key] = d_bin[key] * d_cal[key + '_multiplier'] + \
                        d_cal[key + '_offset']
                if 'time' in d_bin:
                    print("Whoa, there's a field named time in fastDAQ_" +
                          str(i_file))
                d['time'] = (range(d[list(d.keys())[0]].size) -
                             d_cal['pretrigger_samples']) * d_cal['dt']
                if 'bindata' in d_bin:
                    print("Whoa, there's a field named bindata in fastDAQ_" +
                          str(i_file))
                d['bindata'] = d_bin
                if 'caldata' in d_bin:
                    print("Whoa, there's a field named caldata in fastDAQ_" +
                          str(i_file))
                d['caldata'] = d_cal

            except:
                print('Failed to load fastDAQ_' + str(i_file))
                traceback.print_exc()
                break

            if i_file == 0:
                event['fastDAQ'] = d
                event['fastDAQ']['multiboards'] = [d]
            else:
                event['fastDAQ']['multiboards'].append(d)

            event['fastDAQ']['loaded'] = True
            i_file += 1

    if ('slowDAQ' in loadlist) or (neglist and '~slowDAQ' not in loadlist):
        try:
            d = rt(os.path.join(event_dir, 'slowDAQ_0.txt'))
            event['slowDAQ'] = d
            event['slowDAQ']['loaded'] = True
        except:
            print('Failed to load slowDAQ_0.txt')
            traceback.print_exc()

    if ('PMTtraces' in loadlist) or (neglist and '~PMTtraces' not in loadlist):
        try:
            d = rb(os.path.join(event_dir, 'PMTtraces.bin'), max_file_size=max_file_size)
            event['PMTtraces'] = d
            event['PMTtraces']['loaded'] = True
        except:
            print('Failed to load PMTtraces')
            traceback.print_exc()


    if ('event' in loadlist) or (neglist and '~event' not in loadlist):
        try:
            with open(os.path.join(event_dir, 'Event.txt'), 'r') as ev_txt:
                ev_str = next(ev_txt)
            ev_dat = ev_str.split()
            event['event']['run_type'] = np.int32(ev_dat[2])
            event['event']['trigger_main'] = np.int32(ev_dat[3])
            event['event']['trigger_cameras'] = np.int32(ev_dat[4])
            event['event']['trigger_PLC'] = np.int32(ev_dat[5])
            event['event']['trigger_slowDAQ'] = np.int32(ev_dat[6])
            event['event']['timestamp'] = np.float64(ev_dat[7])
            event['event']['mstick'] = np.int64(ev_dat[8])
            event['event']['Pset'] = np.float64(ev_dat[9])
            event['event']['livetime'] = np.float64(ev_dat[10])
            event['event']['loaded'] = True
        except:
            print('Failed to load Event')
            traceback.print_exc()

    return event
