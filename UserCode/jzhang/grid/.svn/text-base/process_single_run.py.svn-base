# python sbc_pmttest_processall.py [run_list]
# if run_list is provided, the runs in the list will be processed; otherwise
# the runs in the script will be processed
# $CVMFS_DIR/code/PICOscripts/process_single_run.pl $OPTION_TEXT $RUN $MYSUBDIR/$DATASERIES $CVMFS_DIR/code $MYSUBDIR/output $MYSUBDIR/log > $MYSUBDIR/log/$RUN/process_single_run.log 2>&1

import numpy as np
import SBCcode as sbc
import os
import sys
import argparse


def parse_process_opt(argv):
    parser = argparse.ArgumentParser(description='Process SBC data.')

    # data selections
    parser.add_argument("--run", nargs=1, help="runs to process")
    # processing options
    parser.add_argument("--event", action='store_true')
    parser.add_argument("--pmtfda", action='store_true')
    parser.add_argument("--pmt", action='store_true')
    parser.add_argument("--dytran", action='store_true')
    parser.add_argument("--acoustic", action='store_true')
    parser.add_argument("--history", action='store_true')
    parser.add_argument("--timing", action='store_true')
    parser.add_argument("--images", action='store_true')
    parser.add_argument("--pmtpf", action='store_true')

    # post processing
    parser.add_argument("--tar_runs")
    parser.add_argument("--verbose")

    # positional arguments
    parser.add_argument("data_dir")
    parser.add_argument("code_dir")
    parser.add_argument("output_dir")
    parser.add_argument("log_dir")

    args = parser.parse_args(argv)
    # print(args)

    runlist = []
    if args.run:
        runlist = args.run

    # process_list=['event', 'pmtfda', 'pmt', 'dytran', 'acoustic', 'history', 'timing', 'images', 'pmtpf']
    process_list = []
    if args.event:
        process_list.append('event')
    if args.pmtfda:
        process_list.append('pmtfda')
    if args.pmt:
        process_list.append('pmt')
    if args.dytran:
        process_list.append('dytran')
    if args.acoustic:
        process_list.append('acoustic')
    if args.history:
        process_list.append('history')
    if args.timing:
        process_list.append('timing')
    if args.images:
        process_list.append('images')
    if args.pmtpf:
        process_list.append('pmtpf')

    data_dir = args.data_dir
    code_dir = args.code_dir
    output_dir = args.output_dir
    log_dir = args.log_dir

    return runlist, process_list, data_dir, code_dir, output_dir, log_dir


def main(argv):
    runlist, process_list, data_dir, code_dir, output_dir, log_dir = parse_process_opt(argv)

    for runname in runlist:
        runid_str = runname.split('_')
        runid = np.int32(runid_str)
        sbc.psr2(os.path.join(data_dir, runname),
                 dataset='SBC-2017',
                 recondir=os.path.join(output_dir, runname),
                 process_list=process_list)
        # process_list=['event', 'pmtfda', 'pmt', 'dytran', 'acoustic', 'history', 'timing', 'images'])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
