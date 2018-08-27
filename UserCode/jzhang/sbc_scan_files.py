from __future__ import print_function
import os
import re
import numpy as np
import sys
import argparse
from operator import itemgetter, attrgetter, methodcaller


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--new", action='store_true', help='new runs starting from a run')
    parser.add_argument("--missing", action='store_true', help='scan missing files')
    parser.add_argument("--all", action='store_true', help='list all runs')
    parser.add_argument("--run_type", action='store_true', help='list run types')
    parser.add_argument("--process", action='store_true', help='list processing info')

    parser.add_argument("start_run", nargs='?')
    parser.add_argument("end_run", nargs='?')

    args = parser.parse_args()
    # print(args)

    start_run = '20170619_3'
    end_run = '20991212_99'
    data_dir = '/bluearc/storage/SBC-17-data'

    if args.start_run:
        start_run = args.start_run
    if args.end_run:
        end_run = args.end_run

    if args.new:
        runlist = get_runlist(data_dir, start_run, end_run)
        for run in runlist:
            print(run)

    if args.missing:
        scan_missing_files(data_dir, start_run)

    if args.all:
        runlist = get_runlist(data_dir, '20170101_0')
        for run in runlist:
            print(run)

    if args.run_type:
        scan_PTR(data_dir, start_run)

    if args.process:
        process_status()


def run_to_sn(runs):
    if type(runs) == str:
        runlist = [runs]
    elif type(runs) == list:
        runlist = runs
    runid = [np.int64(run.split('_')) for run in runlist]
    runid = np.int64(runid)

    runsn = runid[:, 0] * 1000 + runid[:, 1]

    return runsn


def get_runlist(data_dir, start_run='20170101_0', end_run='20991212_99'):
    runlist = os.listdir(data_dir)
    runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None)
                                and os.path.isdir(os.path.join(data_dir, fn))
                                and (len(os.listdir(os.path.join(data_dir, fn))) > 0), runlist)
    runlist = list(runlist)
    # print(runlist)

    sn0 = run_to_sn(start_run)
    sn1 = run_to_sn(end_run)
    runsn = run_to_sn(runlist)

    # sort
    ix = np.argsort(runsn)
    runsn = runsn[ix]
    runlist = [runlist[x] for x in ix]

    ix = np.nonzero(np.logical_and(sn0 <= runsn, runsn <= sn1))[0]

    return [runlist[x] for x in ix]


def scan_missing_files(data_dir, start_run='20170101_0', end_run='20991212_99'):
    run_files = ['DAQ30l_Setup.xml', 'DAQversion.txt', 'RunParameters.txt']
    evt_files = ['Event.txt', 'PLClog.txt', 'PMTtraces.bin', 'fastDAQ_0.bin', 'fastDAQ_0_cal.txt', 'slowDAQ_0.txt']

    runlist = get_runlist(data_dir, start_run, end_run)
    runsn = run_to_sn(runlist)
    t0 = run_to_sn(start_run)[0]
    for irun in range(0, len(runlist)):
        if runsn[irun] >= t0:

            runname = runlist[irun]
            this_run_files = run_files + [runname + '.txt']
            file_list = os.listdir(os.path.join(data_dir, runname))

            if len(file_list) < 1:
                print("Empty run " + os.path.join(data_dir, runname))
                continue

            evt_list = [ev for ev in file_list
                        if (re.search(r"^\d+$", ev)
                            and os.path.isdir(os.path.join(data_dir, runname, ev)))]
            ix = np.argsort(np.int32(evt_list))
            evt_list = [evt_list[x] for x in ix]

            other_list = [f for f in file_list if os.path.isfile(os.path.join(data_dir, runname, f))]

            for f in this_run_files:
                if not f in set(other_list):
                    print('Missing ' + os.path.join(data_dir, runname, f))

            for ev in evt_list:
                files = os.listdir(os.path.join(data_dir, runname, ev))
                for f in evt_files:
                    if not f in set(files):
                        print('Missing ' + os.path.join(data_dir, runname, ev, f))


# Scan pressure, temperature, and run type
def scan_PTR(data_dir, start_run='20170101_0', end_run='20991212_99'):
    runlist = get_runlist(data_dir, start_run, end_run)
    runsn = run_to_sn(runlist)
    t0 = run_to_sn(start_run)[0]
    for irun in range(0, len(runlist)):
        if runsn[irun] >= t0:

            runname = runlist[irun]
            file_list = os.listdir(os.path.join(data_dir, runname))

            if len(file_list) < 1:
                # print("Empty run " + os.path.join(data_dir, runname))
                continue

            evt_list = [ev for ev in file_list
                        if (re.search(r"^\d+$", ev)
                            and os.path.isdir(os.path.join(data_dir, runname, ev)))]

            filename = os.path.join(data_dir, runname, 'DAQ30l_Setup.xml')
            with open(filename, 'r') as f:
                read_data = f.readlines()
                for line in read_data:
                    if re.search(r'<Run_type>(.*)</Run_type>', line):
                        m = re.search(r'<Run_type>(.*)</Run_type>', line)
                        print(runname + '  ' + str(len(evt_list)) + '  ' + m.group(1))


def process_status():
    run_list_dir = '/pnfs/coupp/persistent/runlists/SBC-17'
    skip_pmt_lists = ['skippmt_1500653829', 'skippmt_1502827534']
    pmt_lists = ['pmt_1500770663', 'pmt_1502755802']

    no_pmt_runs = []
    for f in skip_pmt_lists:
        fh = open(os.path.join(run_list_dir, f), 'r')
        no_pmt_runs += fh.readlines()
        fh.close()

    pmt_runs = []
    for f in pmt_lists:
        fh = open(os.path.join(run_list_dir, f), 'r')
        pmt_runs += fh.readlines()
        fh.close()

    no_pmt_runs = set(no_pmt_runs)
    pmt_runs = set(pmt_runs).difference(no_pmt_runs)

    status = ([(run_to_sn(run.strip()), run.strip(), 'no pmt') for run in no_pmt_runs] +
              [(run_to_sn(run.strip()), run.strip(), '') for run in pmt_runs])
    # print(status)
    status = sorted(status, key=itemgetter(0))
    for item in status:
        print(item[1], item[2])


if __name__ == '__main__':
    sys.exit(main())
