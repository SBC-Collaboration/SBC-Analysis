#!/usr/bin/env bash

echo Start `date`
echo Site:${GLIDEIN_ResourceName}
echo "The worker node is " `hostname` " OS: " `uname -a`
echo "The user id is " `whoami`
echo "The output of id is " `id`

echo
echo "df -h"
df -h

# Request for addition to PRODUCTS path - RITM0074989

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
setup ifdhc
setup sam_web_client

echo
echo "source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setups"
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setups
setup root v5_34_12 -q nu:e4:prof
#setup root v6_08_06g -q nu:e14:prof
echo "Setting up gcc v6_3_0"
setup gcc v6_3_0
echo "setup python v2_7_13d"
setup python v2_7_13d


prev_time=`date +%s`
t0=$prev_time

while [[ $1 =~ -- ]]
do
  OPTION_TEXT="$OPTION_TEXT $1"
  shift
done


CVMFS_DIR="/cvmfs/coupp.opensciencegrid.org"
RUNLIST_DIR="/pnfs/coupp/persistent/runlists"

EXE_DIR=$1 # location of SBCcode.tar
DATA_DIR=$2 # eg. /pnfs/coupp/storage/SBC-17-data
CODE_DIR=$3 # CODE_DIR="$CVMFS_DIR/code"
OUTPUT_DIR=$4 # OUTPUT_DIR="/pnfs/coupp/persistent/grid_output"
LOG_DIR=$5 # not used since OUTPUT_DIR is at a higher level
DATASERIES_SHORT=$6
RUN_LIST_FILE=$7 # full path
SKIP_PMT=${8// } # trim all spaces

DATASERIES=${DATA_DIR##*/} # eg. SBC-17-data


if [ -d $_CONDOR_SCRATCH_DIR ]; then
	MYSUBDIR=$_CONDOR_SCRATCH_DIR/process_single_run_$$
else
	MYSUBDIR=$OSG_WN_TMP/process_single_run_$$
fi
echo "MYSUBDIR: $MYSUBDIR"
mkdir $MYSUBDIR



# copy SBCcode
echo
echo "Copying SBCcode"
if ifdh ls $EXE_DIR/SBCcode.tar; then
    ifdh cp -D $EXE_DIR/SBCcode.tar $MYSUBDIR
    cd $MYSUBDIR
    tar xf SBCcode.tar
    echo "ls -l $MYSUBDIR"
    ls -l $MYSUBDIR
else
    echo "$EXE_DIR/SBCcode.tar not found, exit"
    exit 1
fi


cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo SBCcode COPY IN TOOK $minutes MINUTES
prev_time=$cur_time



# export PATH="$CODE_DIR/sbc_python2.7/bin:$PATH"
echo
echo "python path is `which python`"
export PYTHONPATH=$MYSUBDIR:$CODE_DIR/sbc_python2.7.13d/lib/python2.7/site-packages:$PYTHONPATH
echo "PYTHONPATH is $PYTHONPATH"

# Get runname and setup directories
echo
echo "Going to look in $RUN_LIST_FILE for run_list"
TEMP_RUN_LIST_FILE=$MYSUBDIR/$(basename $RUN_LIST_FILE)
echo "ifdh cp $RUN_LIST_FILE $TEMP_RUN_LIST_FILE"
ifdh cp $RUN_LIST_FILE $TEMP_RUN_LIST_FILE
echo "ls -l $MYSUBDIR"
ls -l $MYSUBDIR
LINE=`expr $PROCESS + 1`
RUN=`head -$LINE $TEMP_RUN_LIST_FILE | tail -1`
echo RUN is $RUN



mkdir $MYSUBDIR/$DATASERIES
mkdir $MYSUBDIR/output
mkdir $MYSUBDIR/output/$RUN
mkdir $MYSUBDIR/log
mkdir $MYSUBDIR/log/$RUN
mkdir $MYSUBDIR/tmp
cd $MYSUBDIR


cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo PREAMBLE TOOK $minutes MINUTES
prev_time=$cur_time



# copy data
echo; echo;
if ifdh ls $DATA_DIR/$RUN.tgz; then
    mkdir $MYSUBDIR/$DATASERIES/$RUN
    echo "ifdh cp $DATA_DIR/$RUN.tgz $MYSUBDIR/$DATASERIES/$RUN/$RUN.tgz"
    ifdh cp $DATA_DIR/$RUN.tgz $MYSUBDIR/$DATASERIES/$RUN/$RUN.tgz
    ls $MYSUBDIR/$DATASERIES/$RUN
    tar -xzf $MYSUBDIR/$DATASERIES/$RUN/$RUN.tgz -C $MYSUBDIR/$DATASERIES/$RUN
    if [ -e $MYSUBDIR/$DATASERIES/$RUN/$RUN ]; then
        echo "Found run within run, moving it out"
        chmod -R +w $MYSUBDIR
        mv $MYSUBDIR/$DATASERIES/$RUN $MYSUBDIR/$DATASERIES/$RUN\_tmp
        mv $MYSUBDIR/$DATASERIES/$RUN\_tmp/$RUN $MYSUBDIR/$DATASERIES
        rm -r $MYSUBDIR/$DATASERIES/$RUN\_tmp
    fi
else
    if ifdh ls $DATA_DIR/$RUN.tar.bz2; then
        mkdir $MYSUBDIR/$DATASERIES/$RUN
        echo "ifdh cp $DATA_DIR/$RUN.tar.bz2 $MYSUBDIR/$DATASERIES/$RUN/$RUN.tar.bz2"
        ifdh cp $DATA_DIR/$RUN.tar.bz2 $MYSUBDIR/$DATASERIES/$RUN/$RUN.tar.bz2
        ls $MYSUBDIR/$DATASERIES/$RUN
        tar -xjf $MYSUBDIR/$DATASERIES/$RUN/$RUN.tar.bz2 -C $MYSUBDIR/$DATASERIES/$RUN
        if [ -e $MYSUBDIR/$DATASERIES/$RUN/$RUN ]; then
            echo "Found run within run, moving it out"
            chmod -R +w $MYSUBDIR
            mv $MYSUBDIR/$DATASERIES/$RUN $MYSUBDIR/$DATASERIES/$RUN\_tmp
            mv $MYSUBDIR/$DATASERIES/$RUN\_tmp/$RUN $MYSUBDIR/$DATASERIES
            rm -r $MYSUBDIR/$DATASERIES/$RUN\_tmp
        fi
    else
    	echo "Run not in $DATA_DIR"
    fi
fi
echo
echo "ls -larth $MYSUBDIR/$DATASERIES/$RUN"
ls -larth $MYSUBDIR/$DATASERIES/$RUN

#echo "Computing md5sum"
#pushd $MYSUBDIR/$DATASERIES/$RUN
#md5sum *.* > $MYSUBDIR/log/$RUN/$RUN.md5
#dirlist=`ls -d */ | sort -V`
#for d in $dirlist; do
#	md5sum $d* >> $MYSUBDIR/log/$RUN/$RUN.md5
#done
#popd
#echo "Done computing md5sum"


cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo DATA COPY IN TOOK $minutes MINUTES
prev_time=$cur_time



echo; echo;
echo LOOKING FOR CODE DIRECTORY
echo "ls -larth $CODE_DIR"
ls -larth $CODE_DIR


cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo CODE LS TOOK $minutes MINUTES
prev_time=$cur_time



echo; echo;
if [ -e $CODE_DIR/Matlab/LookupTables/$DATASERIES_SHORT ]; then
    echo "ls -larth $CODE_DIR/Matlab/LookupTables/$DATASERIES_SHORT"
    ls -larth $CODE_DIR/Matlab/LookupTables/$DATASERIES_SHORT
fi


cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo LOOKUP LS IN TOOK $minutes MINUTES
prev_time=$cur_time


echo; echo;
echo "ls -larth $CVMFS_DIR/MCR"
ls -larth $CVMFS_DIR/MCR


cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo MCR LS  TOOK $minutes MINUTES
prev_time=$cur_time



chmod -R +w $MYSUBDIR

export PERL5LIB=$CODE_DIR:$CODE_DIR/PICOscripts:$PERL5LIB
MCRROOT=$CVMFS_DIR/MCR/current
MCRJRE=$MCRROOT/sys/java/jre/glnxa64/jre/lib/amd64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/runtime/glnxa64:$MCRROOT/bin/glnxa64:$MCRROOT/sys/os/glnxa64:$MCRROOT/sys/opengl/lib/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRJRE/native_threads:$MCRJRE/server:$MCRJRE/client:$MCRJRE
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CODE_DIR/AutoBub2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CODE_DIR/AutoBub3
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CODE_DIR/AutoBub3hs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CODE_DIR/jartracker
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export XAPPLRESDIR=$MCRROOT/X11/app-defaults
export MCR_CACHE_ROOT=$MYSUBDIR/tmp

#printenv

###############################################################################

echo; echo;
echo "Doing analysis"
echo "python $CODE_DIR/SBCcode/UserCode/jzhang/grid/process_single_run.py $OPTION_TEXT $RUN $MYSUBDIR/$DATASERIES $CODE_DIR $MYSUBDIR/output $MYSUBDIR/log > $MYSUBDIR/log/$RUN/process_single_run.log 2>&1"
#python $CODE_DIR/SBCcode/UserCode/jzhang/grid/process_single_run.py $OPTION_TEXT --run $RUN $MYSUBDIR/$DATASERIES $CODE_DIR $MYSUBDIR/output $MYSUBDIR/log > $MYSUBDIR/log/$RUN/process_single_run.log 2>&1

#if [[ -z $SKIP_PMT ]]; then
#    OPTION_TEXT="--event --acoustic --dytran --history --images --pmtfda --pmt --timing"
#else
#    OPTION_TEXT="--event --acoustic --dytran --history --images"
#fi
python $CODE_DIR/SBCcode/UserCode/jzhang/grid/process_single_run.py \
    $OPTION_TEXT \
    --run $RUN $MYSUBDIR/$DATASERIES $CODE_DIR $MYSUBDIR/output $MYSUBDIR/log > $MYSUBDIR/log/$RUN/process_single_run.log 2>&1

#if [ -e $MYSUBDIR/log/$RUN/processing_failed ]
#then
#    $CODE_DIR/PICOscripts/process_single_run.pl $OPTION_TEXT $RUN $MYSUBDIR/$DATASERIES $CODE_DIR $MYSUBDIR/output $MYSUBDIR/log >> $MYSUBDIR/log/$RUN/process_single_run.log 2>&1
#fi

echo; echo;
cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo ANALYSIS TOOK $minutes MINUTES
prev_time=$cur_time



echo; echo;
echo "CHECKING DIRECTORY AT $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN:"
if ! ifdh ls $OUTPUT_DIR/$DATASERIES_SHORT &> /dev/null; then
    ifdh mkdir $OUTPUT_DIR/$DATASERIES_SHORT
    echo "Making $OUTPUT_DIR/$DATASERIES_SHORT directory!"
fi
if ! ifdh ls $OUTPUT_DIR/$DATASERIES_SHORT/output &> /dev/null; then
    ifdh mkdir $OUTPUT_DIR/$DATASERIES_SHORT/output
    echo "Making $OUTPUT_DIR/$DATASERIES_SHORT/output directory!"
fi
if ! ifdh ls $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN &> /dev/null; then
    ifdh mkdir $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN
    echo "Making $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN directory!"
fi

OUTFILES=$(ls $MYSUBDIR/output/$RUN/*.bin)
for FF in $OUTFILES; do
    ifdh rm $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN/$(basename $FF) &> /dev/null
    echo COPYING $FF
    ifdh cp -D $FF $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN/
done

#if ifdh ls $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN &> /dev/null; then
#    echo "Removing $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN"
#    ifdh rm $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN/* 2>&1
#    ifdh rmdir $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN 2>&1
#fi
#
#echo COPYING $MYSUBDIR/output/$RUN $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN
##chmod -R +w $OUTPUT_DIR
#ifdh cp -r $MYSUBDIR/output/$RUN $OUTPUT_DIR/$DATASERIES_SHORT/output/$RUN 2>&1
##chmod -R +w $OUTPUT_DIR/$RUN

cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo OUTPUT COPY OUT TOOK $minutes MINUTES
prev_time=$cur_time



echo; echo;
echo "CHECKING DIRECTORY AT $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN:"
if ! ifdh ls $OUTPUT_DIR/$DATASERIES_SHORT/log &> /dev/null; then
    ifdh mkdir $OUTPUT_DIR/$DATASERIES_SHORT/log
    echo "Making $OUTPUT_DIR/$DATASERIES_SHORT/log directory!"
fi
if ! ifdh ls $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN &> /dev/null; then
    ifdh mkdir $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN
    echo "Making $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN"
fi
LOGFILES=$(ls $MYSUBDIR/log/$RUN/*.log)
for FF in $LOGFILES; do
    ifdh rm $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN/$(basename $FF) &> /dev/null
    echo "Copying $FF"
    ifdh cp -D $FF $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN
done


#if ifdh ls $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN &> /dev/null; then
#    echo "Removing $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN"
#    ifdh rm $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN/* 2>&1
#    ifdh rmdir $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN 2>&1
#fi
#echo "ifdh cp -r $MYSUBDIR/log/$RUN $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN 2>&1"
#ifdh cp -r $MYSUBDIR/log/$RUN $OUTPUT_DIR/$DATASERIES_SHORT/log/$RUN 2>&1


cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo LOG COPY OUT TOOK $minutes MINUTES
prev_time=$cur_time


echo; echo;
echo "Directories and files on the processing node"
du -h --max-depth=1 $MYSUBDIR
echo
echo $MYSUBDIR:
ls -larth $MYSUBDIR/
echo
echo $MYSUBDIR/output/$RUN:
ls -larth $MYSUBDIR/output/$RUN
echo
echo $MYSUBDIR/log/$RUN:
ls -larth $MYSUBDIR/log/$RUN



echo; echo;
echo "Cleaning up"
echo "Changing permissions and removing $MYSUBDIR"
chmod -R a+w $MYSUBDIR 2>&1
rm -fr  $MYSUBDIR/* 2>&1



echo; echo;
cur_time=`date +%s`
minutes=`echo "scale=2; ($cur_time-$prev_time)/60.0" | bc`
echo POSTAMBLE TOOK $minutes MINUTES
minutes=`echo "scale=2; ($cur_time-$t0)/60.0" | bc`
echo ENTIRE PROCESSING TOOK $minutes MINUTES
