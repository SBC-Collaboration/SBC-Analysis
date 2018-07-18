#!/usr/bin/env bash

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
setup ifdhc
setup sam_web_client

echo "source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setups.for.cvmfs"
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setups.for.cvmfs
setup root v5_34_12 -q nu:e4:prof
#setup root v6_08_06g -q nu:e14:prof
echo "Setting up gcc v6_3_0"
setup gcc v6_3_0
echo "setup python v2_7_13d"
setup python v2_7_13d

prev_time=`date +%s`
t0=$prev_time


CVMFS_DIR="/cvmfs/coupp.opensciencegrid.org"
EXE_DIR="/cvmfs/coupp.opensciencegrid.org/code" # not used by the wrapper
# no trailing slash
DATA_DIR="/bluearc/storage/SBC-17-data"
CODE_DIR="/bluearc/storage/recon/devel/code"
# CODE_DIR="/coupp/data/storage/recon/devel/code" # need to publish code or put on pnfs
OUTPUT_DIR="/bluearc/storage/recon/devel/SBC-17/output"
LOG_DIR="/bluearc/storage/recon/devel/SBC-17/log"

PNFS_OUTPUT_DIR="/pnfs/coupp/persistent/grid_output"
PNFS_LOG_DIR="/pnfs/coupp/persistent/grid_output"


#RUN_LIST_FILE="/bluearc/storage/SBC-17-tgz/SBC_proc_runlist"
RUN_LIST_FILE="/bluearc/storage/SBC-17-tgz/listlocal_1499715541"

DATASERIES=${DATA_DIR##*/} # eg. SBC-17-data
DATASERIES_SHORT="SBC-17"

echo "python path is `which python`"
export PYTHONPATH=$CODE_DIR:$CODE_DIR/sbc_python2.7.13d/lib/python2.7/site-packages:$PYTHONPATH
echo "PYTHONPATH is $PYTHONPATH"

proc_runs=$(cat $RUN_LIST_FILE)
for RUN in $proc_runs; do
    echo "Processing $RUN"
    if [ ! -d $LOG_DIR/$RUN ]; then
        mkdir $LOG_DIR/$RUN
    else 
        rm -rf $LOG_DIR/$RUN
        mkdir $LOG_DIR/$RUN
    fi
    
    if [ ! -d $OUTPUT_DIR/$RUN ]; then
        mkdir $OUTPUT_DIR/$RUN
    else 
        rm -rf $OUTPUT_DIR/$RUN
        mkdir $OUTPUT_DIR/$RUN
    fi
#    python $CODE_DIR/SBCcode/Tools/FermiGrid/process_single_run.py \
#    --event --pmtfda --pmt --dytran --acoustic --history --timing --images \
#    --run $RUN $DATA_DIR $CODE_DIR $OUTPUT_DIR $LOG_DIR > $LOG_DIR/$RUN/process_single_run.log 2>&1
    python $CODE_DIR/SBCcode/Tools/FermiGrid/process_single_run.py \
    --event --acoustic --dytran --history --images \
    --run $RUN $DATA_DIR $CODE_DIR $OUTPUT_DIR $LOG_DIR > $LOG_DIR/$RUN/process_single_run1.log 2>&1
        python $CODE_DIR/SBCcode/Tools/FermiGrid/process_single_run.py \
    --pmtfda --pmt --acoustic --timing \
    --run $RUN $DATA_DIR $CODE_DIR $OUTPUT_DIR $LOG_DIR > $LOG_DIR/$RUN/process_single_run2.log 2>&1
done

#for RUN in $proc_runs; do
#    echo "CHECKING DIRECTORY AT $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN:"
#    if ! ifdh ls $PNFS_OUTPUT_DIR/$DATASERIES_SHORT &> /dev/null; then
#        ifdh mkdir $PNFS_OUTPUT_DIR/$DATASERIES_SHORT
#        echo "Making $PNFS_OUTPUT_DIR/$DATASERIES_SHORT directory!"
#    fi
#    if ! ifdh ls $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output &> /dev/null; then
#        ifdh mkdir $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output
#        echo "Making $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output directory!"
#    fi
#    if ifdh ls $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN &> /dev/null; then
#        echo "Removing $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN"
#        ifdh rm $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN/* 2>&1
#        ifdh rmdir $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN 2>&1
#    fi
#
#
#    echo COPYING  $OUTPUT_DIR/$RUN $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN
#    #chmod -R +w $OUTPUT_DIR
#    ifdh cp -r $OUTPUT_DIR/$RUN $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN 2>&1
#    #chmod -R +w $OUTPUT_DIR/$RUN
#
#    echo "CHECKING DIRECTORY AT $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN:"
#    if ! ifdh ls $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log &> /dev/null; then
#        ifdh mkdir $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log
#        echo "Making $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log directory!"
#    fi
#    if ifdh ls $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN &> /dev/null; then
#        echo "Removing $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN"
#        ifdh rm $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN/* 2>&1
#        ifdh rmdir $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN 2>&1
#    fi
#    echo "ifdh cp -r $LOG_DIR/$RUN $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN 2>&1"
#    ifdh cp -r $LOG_DIR/$RUN $PNFS_OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN 2>&1
#
#    cur_time=`date +%s`
#    minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
#    echo POSTAMBLE TOOK $minutes MINUTES
#    minutes=`echo "scale=2; ($cur_time-$t0)/60.0" | bc`
#    echo ENTIRE PROCESSING TOOK $minutes MINUTES
#done
