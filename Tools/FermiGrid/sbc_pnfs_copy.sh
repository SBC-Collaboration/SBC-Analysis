#!/usr/bin/env bash
# cron job file to tar data to pnfs

exec 11>/home/coupp/syncSBCpnfs
if ! flock -n 11; then
	# echo "another instance is running";
	exit 1
else
    exec 9>/home/coupp/syncSBC
	if ! flock -n 9; then
	    # echo "Data syncing is running"
	    exit 2
	else
	    flock -u 9
	fi
fi


data_dir="/bluearc/storage/SBC-17-data"
tgz_dir="/bluearc/storage/SBC-17-tgz"
pnfs_dir="/pnfs/coupp/storage/SBC-17-data"
output_dir="/pnfs/coupp/persistent/grid_output/SBC-17/output"

# all runs in data_dir
if cd $data_dir; then
    data_runs=`ls -d 2017*_*/ | sed 's#/##' | sort -V`
else
    exit 1
fi

# all runs in pnfs_dir
if cd $pnfs_dir; then
    pnfs_runs=`ls 2017*_*.tgz | sed 's#.tgz##' | sort -V`
else
    exit 1
fi

#- for processing
if cd $output_dir; then
    out_runs=`ls -d 2017*_*/ | sed 's#/##' | sort -V`
else
    exit 1
fi
#-

# runs that are NOT uptodate as reported by rsync
if [ -f $tgz_dir/rsync.log ]; then
    rm -rf $tgz_dir/rsync.log
fi
rsync -artv --stats --chmod=a-w,a+r,Da+x  --dry-run \
    pico@dm.phys.northwestern.edu:/mnt/XENON_DAQ/SBC-17-data \
    /bluearc/storage/ > $tgz_dir/rsync.log
if [ $? -eq 0 ]; then
    new_runs=`cat $tgz_dir/rsync.log | grep '^SBC' | sed 's#.*data/\([0-9]*_[0-9]*\).*#\1#' | uniq`
else
    exit 1
fi
# new_runs may contain incomplete runs, so remove them from processing or taring

#- for processing
proc_runs=${data_runs#*20170619_2} # remove earlier runs
for run in $new_runs; do
    proc_runs=`echo $proc_runs | sed "s#$run##"`
done
for run in $out_runs; do
    proc_runs=`echo $proc_runs | sed "s#$run##"`
done
proc_runs=`echo $proc_runs | xargs` # trim spaces
echo -n > $tgz_dir/SBC_proc_runlist
if [[ ! -z $proc_runs ]]; then
    for run in $proc_runs; do
        echo $run >> $tgz_dir/SBC_proc_runlist
    done
fi
#-

# for taring
for run in $new_runs; do
    data_runs=`echo $data_runs | sed "s#$run##"`
done

for run in $pnfs_runs; do
    data_runs=`echo $data_runs | sed "s#$run##"`
done

data_runs=`echo $data_runs | xargs` # trim spaces


if [[ ! -z $data_runs ]]; then
    echo `date`
    echo "Runs to tar:"
    echo $data_runs
    for run in $data_runs; do
      if [ ! -f $tgz_dir/$run.tgz ]; then
#        echo "$data_dir/$run/$run.tgz exists, removing"
#        chmod +w $data_dir/$run
#        rm -rf $data_dir/$run/$run.tgz
        cd $data_dir/$run
        echo "Taring $tgz_dir/$run.tgz"
        # touch $data_dir/$run/$run.tgz
        tar --ignore-failed-read --exclude=*.tgz -czf $tgz_dir/$run.tgz *
        # chmod -w $data_dir/$run
      fi

#    if [ -f $pnfs_dir/$run.tgz ]; then
#        echo "$data_dir/$run/$run.tgz exists, removing"
#        chmod +w $pnfs_dir/$run.tgz
#        rm -rf $pnfs_dir/$run.tgz
#    fi

    echo "ifdh cp -D $tgz_dir/$run.tgz $pnfs_dir"
    ifdh cp -D $tgz_dir/$run.tgz $pnfs_dir
    if cmp $tgz_dir/$run.tgz $pnfs_dir/$run.tgz &> /dev/null  # Suppress output.
    then
        echo "$run.tgz copied and verified"
        rm -rf $tgz_dir/$run.tgz
    else
        echo "$run.tgz failed"
    fi
    done
fi
#echo $data_runs
#echo ${data_runs%2017*_*} # remove last run
#echo $pnfs_runs

# cat sync.log | grep '^SBC' | grep -v 'uptodate' | sed 's#.*data/\([0-9]*_[0-9]*\).*#\1#' | uniq
