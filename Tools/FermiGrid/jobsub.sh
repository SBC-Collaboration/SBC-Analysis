#!/usr/bin/env bash
#
#cd /bluearc/storage/recon/devel/code/SBCcode/AnalysisModules
#svn up
#cd /bluearc/storage/recon/devel/code/SBCcode/Tools/FermiGrid
#./tar_sbc_code.sh
### OLD RUNS
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_20170705" # where to get the run for the processing script
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499435272" # one run test
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499442343" # 110 runs
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499653396" # 63 runs
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499659072" # re-sub
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499693254" # 13 runs
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499708559" # 24 runs, corrected dytran/acoustic flip
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499757474" # 8 new runs
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499761912" # 27 runs re-sub
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499793466"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499881227"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499911668"
# RUN_LIST_F" # all runs for Piezo2
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1500323827"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1500577322"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1500592889" # fix acoustic dimensions
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1500627181" # all runs for corrected history analysis
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/skippmt_1500653829"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/skippmt_1502827534"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/image_1502941012"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/pmtpulse_1503095464"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/pmtpulse_1503112494" # --expected-lifetime=long
# SKIP_PMT="Y"
# OPTION_TEXT="--event --acoustic --dytran --history --images"
# OPTION_TEXT="--pmtpf"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/pmt_1500770663" # with T0finder2
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/pmt_1502755802"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1507619130"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1509336953"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1509349070"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_SepOct"
# RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1510013290"
# SKIP_PMT=""
# OPTION_TEXT="--event --acoustic --dytran --history --images --pmtfda --pmt --timing" # --expected-lifetime=short
# OPTION_TEXT="--event --acoustic --dytran --history"

# EXE_DIR="/cvmfs/coupp.opensciencegrid.org/code" # not used by the wrapper
EXE_DIR="/pnfs/coupp/persistent/runlists/SBC-17" # store SBCcode.tar
# no trailing slash
DATA_DIR="/pnfs/coupp/storage/SBC-17-data"
CODE_DIR="/cvmfs/coupp.opensciencegrid.org"
# CODE_DIR="/coupp/data/storage/recon/devel/code" # need to publish code or put on pnfs
OUTPUT_DIR="/pnfs/coupp/persistent/grid_output/"
LOG_DIR="/pnfs/coupp/persistent/grid_output"
DATASERIES_SHORT="SBC-17-T0Test3"

RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_SepOct"
SKIP_PMT=""
OPTION_TEXT="--acoustic"
#OPTION_TEXT="--acoustic"
N=`cat $RUN_LIST_FILE | wc -l`
jobsub_submit --group=coupp --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC \
    --OS=SL6 \
    --memory=8000MB \
    --disk=80GB \
    --expected-lifetime=medium \
    -N $N \
    file:///nashome/j/jgresl/Projects/Event-Viewer/SBCcode/Tools/FermiGrid/sbc_grid_wrapper.sh \
    $OPTION_TEXT \
    $EXE_DIR \
    $DATA_DIR \
    $CODE_DIR \
    $OUTPUT_DIR \
    $LOG_DIR \
    $DATASERIES_SHORT \
    $RUN_LIST_FILE \
    $SKIP_PMT |& tee jobsub_$$.log