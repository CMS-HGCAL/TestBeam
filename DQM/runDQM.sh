#!/bin/bash
runType=$1
runNumber=$2

dqmOutDir=/tmp/dqmOutDir
dqmPlotDirBase=/var/www/html/dqm

outputFile=${runType}_Output_`printf %06d ${runNumber}`_Reco.root

file=${dqmOutDir}/${outputFile}
nSpills=15

dqmPlotDir=${dqmPlotDirBase}/${runNumber}
mkdir -p $dqmPlotDir/Detailed || exit 


cmsRun test_cfg.py runType=${runType} chainSequence=4  outputFolder=${dqmOutDir} runNumber=${runNumber} maxEvents=-1 > /dev/null


#./DQM/dumpPlots.sh ${dqmOutDir}/${outputFile} ${runNumber}

root -l -b -q "DQM/DumpPlotsReco.C++(\"$file\", \"$dqmPlotDir\", $runNumber, $nSpills)"





