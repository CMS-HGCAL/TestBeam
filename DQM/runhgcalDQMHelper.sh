#!/bin/bash

# echo "1: $1 2: $2"

# CWD=$(pwd)
# cd $1
# SCRAM_ARCH=slc6_amd64_gcc530
# eval `scram runtime -sh`
# # cmsRun test_cfg.py print dataFolder=$2 outputFolder=$3 commonPrefix=$4 runNumber=$5 runType=$6 chainSequence=$7 pedestalsHighGain=$8 pedestalsLowGain=$9
# COMMAND_TO_RUN= "cmsRun test_cfg.py $2"
# echo "Command to run: "
# echo $COMMAND_TO_RUN
# $COMMAND_TO_RUN
# cd $CWD

export SCRAM_ARCH=slc6_amd64_gcc530
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $1 && cd .. && eval `scram runtime -sh` && cd $1
python runhgcalDQM.py
