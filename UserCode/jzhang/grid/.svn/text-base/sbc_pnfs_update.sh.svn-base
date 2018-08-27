#!/usr/bin/env bash
data_dir="/bluearc/storage/SBC-17-data"
cd $data_dir
pnfs_dir="/pnfs/coupp/storage/SBC-17-data"

#RUN_LIST="/pnfs/coupp/persistent/runlists/SBC-17/list_1499442343"
#for LINE in 19 20 30 34 35 37 43 109
#do
#	LINE=`expr $LINE + 1`
#    run=`head -$LINE $RUN_LIST | tail -1`
#    echo $run
#done


for run in \
20170623_3 \
20170624_4 \
20170625_0 \
20170627_0 \
20170629_0 \
20170630_1 \
20170630_8 \
20170702_4 \
20170703_7 \
20170704_3 \
20170704_4
do
    if [ -f $data_dir/$run/$run.tgz ]; then
        echo "$data_dir/$run/$run.tgz exists, removing"
        chmod +w $data_dir/$run
        rm -rf $data_dir/$run/$run.tgz
        cd $data_dir/$run
        echo "Retaring $data_dir/$run/$run.tgz"
        touch $data_dir/$run/$run.tgz
        tar --ignore-failed-read --exclude=*.tgz -czf $run.tgz *
        chmod -w $data_dir/$run
    fi

    if [ -f $pnfs_dir/$run.tgz ]; then
        echo "$data_dir/$run/$run.tgz exists, removing"
        chmod +w $pnfs_dir/$run.tgz
        rm -rf $pnfs_dir/$run.tgz
    fi

    echo "ifdh cp -D $data_dir/$run/$run.tgz $pnfs_dir"
    ifdh cp -D $data_dir/$run/$run.tgz $pnfs_dir
    if cmp $data_dir/$run/$run.tgz $pnfs_dir/$run.tgz &> /dev/null  # Suppress output.
    then echo "$run.tgz copied"
    else echo "$run.tgz failed"
    fi
done
