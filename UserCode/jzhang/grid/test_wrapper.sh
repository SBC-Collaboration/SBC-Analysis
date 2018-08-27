#!/usr/bin/env bash
# $msg = `jobsub_submit $args file://$execdir/process_single_run.grid_wrapper $option_text $execdir $remote_data_dir $remote_code_dir $remote_output_dir $remote_log_dir $data_series run_list_$$ 2>/dev/null`;

# test wrapper
test_dir="/coupp/data/home/coupp/sbc_grid_test"
export PROCESS=98
export _CONDOR_SCRATCH_DIR="$test_dir"
# export _CONDOR_SCRATCH_DIR="$(pwd)/test"

# directory arguments taken by the wrapper
# EXE_DIR=$1
# DATA_DIR=$2
# CODE_DIR=$3
# OUTPUT_DIR=$4
# LOG_DIR=$5
# DATASERIES_SHORT=$6
# RUN_LIST_FILE=$7

EXE_DIR=$test_dir # not used by the wrapper inside
# no trailing slash
DATA_DIR="/pnfs/coupp/storage/SBC-17-data"
CODE_DIR="/coupp/data/storage/recon/devel/code" # OK for testing
OUTPUT_DIR="/pnfs/coupp/persistent/grid_output"
LOG_DIR="/pnfs/coupp/persistent/grid_output"
DATASERIES_SHORT="SBC-17"
RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_20170705" # where to get the run for the processing script

# wrapper just ignores the options
$(pwd)/sbc_grid_wrapper.sh \
--event \
--pmtfda \
--pmt \
--dytran \
--acoustic \
--history \
--timing \
--images \
--pmtpf \
$EXE_DIR \
$DATA_DIR \
$CODE_DIR \
$OUTPUT_DIR \
$LOG_DIR \
$DATASERIES_SHORT \
$RUN_LIST_FILE
