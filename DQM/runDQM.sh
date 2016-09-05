#!/bin/bash
runType=$1
runNumber=$2

dataFolder=/afs/cern.ch/work/r/rslu/public/HGC_TB_data_Sep2016/
dqmOutDir=/tmp/dqmOutDir
dqmPlotDirBase=/var/www/html/dqm
dqmPlotDir=${dqmPlotDirBase}/${runNumber}
mkdir -p $dqmPlotDir/Detailed || exit

chainIdentifier=Reco
nSpills=2
chainSequence=5
pedestalsHighGain=CondObjects/data/Ped_HighGain_L8.txt
pedestalsLowGain=CondObjects/data/Ped_LowGain_L8.txt

if [ "PED" = $runType ]; then
    chainIdentifier=Digi
    nSpills=3
    chainSequence=1
    pedestalsHighGain=${dqmOutDir}/Ped_HighGain_L8_`printf %06d ${runNumber}`.txt
    pedestalsLowGain=${dqmOutDir}/Ped_LowGain_L8_`printf %06d ${runNumber}`.txt
fi

outputFile=${runType}_Output_`printf %06d ${runNumber}`_${chainIdentifier}.root
fileToDump=${dqmOutDir}/${outputFile}

cmsRun test_cfg.py print dataFolder=${dataFolder} outputFolder=${dqmOutDir} runNumber=${runNumber} runType=${runType} chainSequence=${chainSequence} nSpills=${nSpills} pedestalsHighGain=${pedestalsHighGain} pedestalsLowGain=${pedestalsLowGain} maxEvents=-1 > /dev/null

#./DQM/dumpPlots.sh ${dqmOutDir}/${outputFile} ${runNumber}

root -l -b -q "DQM/DumpPlots${chainIdentifier}.C++(\"${fileToDump}\", \"${dqmPlotDir}\", ${runNumber}, ${nSpills})"
