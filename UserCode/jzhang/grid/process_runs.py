import argparse
import sys
import process_single_runs as ps
# /grid/fermiapp/products/common/prd/jobsub_client/v1_2_3_2/NULL
import jobsub_submit as js


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Process SBC data.')

    # data selections
    parser.add_argument("--run", nargs='*')
    parser.add_argument("--runlist")
    parser.add_argument("--all", action='store_true')
    parser.add_argument("--update", action='store_true')
    parser.add_argument("--skip_processed", action='store_true')

    parser.add_argument("--sam_data_stream")
    parser.add_argument("--sam_definition")
    parser.add_argument("--sam_run")

    # cpu selections
    parser.add_argument("--grid", action='store_true')
    parser.add_argument("--gridtest", action='store_true')
    parser.add_argument("--osg", action='store_true')

    # processings
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

    args = parser.parse_args()
    print(args)

    # data to process
    runlist = []
    if args.run and len(args.run) > 0:
        runlist = [run for run in args.run if re.match(r"\d{8}_\d+", run)]
    elif args.runlist:
        if os.path.isfile(args.runlist):
            runlist = read_runlist(args.runlist)
    else:
        if args.all:
            runlist, mtime = get_runlist(data_dir)
        elif args.update:
            starttime = get_starttime(output_dir)
            runlist, mtime= get_runlist(data_dir)
            runlist = [runlist[x] for x in range(0, len(mtime)) if mtime[x] > starttime]
        elif arts.skip_processed:
            runlist, mtime = get_runlist(data_dir)
            outlist, mtime = get_runlist(output_dir)
            runlist = list(set(runlist).difference(set(outlist)))

    exclude_list = excluded_runs('SBC-2017')
    runlist = list(set(runlist).difference(set(exclude_list)))

    if len(runlist) < 1:
        print('runlist is empty, no data need to process.')
        return

    # where to process
    wrapper_file = ''
    if args.gridtest:
        submit_grid_jobs(wrapper_file, runlist, output_dir, log_dir, data_series)
        wait_for_grid_jobs()
    elif args.grid:
        submit_grid_jobs(wrapper_file, runlist, output_dir, log_dir, data_series)
        wait_for_grid_jobs()
    else:
        sbc_argv = ['--run'] + runlist + process_list + [data_dir, code_dir, output_dir, log_dir]
        ps.main(sbc_argv)

    return


def read_runlist(filename):
    with io.open(filename, 'r') as f:
        lines = f.readlines()
    runlist = []
    if len(lines) > 0:
        for line in lines:
            runlist.append(re.findall(r"\d{8}_\d+", line))

    return runlist


def get_runlist(data_dir):
    runlist = os.listdir(data_dir)
    runlist = filter(lambda fn: (not re.search('^\d+_\d+$', fn) is None) and
                                os.path.isdir(os.path.join(datadir, fn)),
                     runlist)
    # runlist = filter(lambda fn: os.path.exists(os.path.join(datadir,
    #                                                         *[fn, 'DAQversion.txt'])), runlist)
    mtime_list = [os.path.getmtime(os.path.join(data_dir, run)) for run in runlist]

    return runlist, mtime_list


def get_starttime(output_dir):
    with open(os.path.join(output_dir, '.process_runs_previous_time')) as f:
        last_time = f.readline()

    return int(last_time)


def submit_grid_jobs(wrapper_file, runlist, output_dir, log_dir, data_series):
    exec_dir = "/coupp/app/condor-exec/coupp/process_runs_$$"

    # ($data_dir, $script_dir, $output_dir, $log_dir, $data_series )

    args = ['jobsub_submit', '--group=coupp', '--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC',
            '-N', '1', '--OS=SL6', '--expected-lifetime=short', '--memory=2000MB', '--disk=40GB']
    args = args + [wrapperfile] + process_list + [exec_dir, data_dir, code_dir, output_dir, log_dir, data_series]

    sys.argv = args # has to change sys.argv due to a reference in jobsub_submit
    js.main(sys.argv)


def wait_for_grid_jobs():
    return 0


def excluded_runs(data_series):
    excluded_runs = []
    if data_series == 'sbc-2015':
        excluded_runs = []

    return excluded_runs

if __name__ == '__main__':
    main()