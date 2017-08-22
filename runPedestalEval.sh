#!/bin/bash

run=$1

path="/afs/cern.ch/work/r/rchatter/TestBeam_July/CMSSW_8_0_21/src/HGCal"
CFG="${path}/unpack2017_cfg.py"
TOP=${PWD}

echo ${path}
echo ${CFG}

cp ${CFG} ./cfg.py

cd ${path}
eval `scramv1 runtime -sh` 
cd ${TOP}

eosPath="/eos/cms/store/group/dpg_hgcal/tb_hgcal/july2017/July2017_TB_data_orm"
echo "Download ${eosPath}/HexaData_Run${run}.raw"

# THIS IS DO DOWNLOAD THE TIMING TXT FILES WHICH YOU DON'T NEED
 echo "Download ${eosPath}/HexaData_Run${run}_TIMING_RDOUT_ORM0.txt"
 xrdcp -f root://eoscms.cern.ch/${eosPath}/HexaData_Run${run}_TIMING_RDOUT_ORM0.txt HexaData_Run${run}_TIMING_RDOUT_ORM0.txt
 echo "Download ${eosPath}/HexaData_Run${run}_TIMING_RDOUT_ORM1.txt"
 xrdcp -f root://eoscms.cern.ch/${eosPath}/HexaData_Run${run}_TIMING_RDOUT_ORM1.txt HexaData_Run${run}_TIMING_RDOUT_ORM1.txt
 echo "Download ${eosPath}/HexaData_Run${run}_TIMING_RDOUT_ORM2.txt"
 xrdcp -f root://eoscms.cern.ch/${eosPath}/HexaData_Run${run}_TIMING_RDOUT_ORM2.txt HexaData_Run${run}_TIMING_RDOUT_ORM2.txt

cmsRun unpack2017_cfg.py runNumber=${run} dataFolder="./" outputFolder="./"

# THIS IS TO SAVE OUTPUT ON EOS : PLEASE DON'T DO THAT
# eosOutputPath="/eos/cms/store/group/dpg_hgcal/tb_hgcal/july2017/pedestalFiles"
# 
# echo "Upload ${eosOutputPath}/PedestalOutput_${run}.root"
# xrdcp -f PedestalOutput_${run}.root root://eoscms.cern.ch/${eosOutputPath}/PedestalOutput_${run}.root
# echo "Upload ${eosOutputPath}/pedestalHG_${run}.txt"
# xrdcp -f pedestalHG_${run}.txt root://eoscms.cern.ch/${eosOutputPath}/pedestalHG_${run}.txt
# echo "Upload ${eosOutputPath}/pedestalLG_${run}.txt"
# xrdcp -f pedestalLG_${run}.txt root://eoscms.cern.ch/${eosOutputPath}/pedestalLG_${run}.txt
# 
# eosOutputPath="/eos/cms/store/group/dpg_hgcal/tb_hgcal/july2017/HGCalTBSkiroc2CMS"
# echo "Upload ${eosOutputPath}/cmsswEvents_Run${run}.root"
# xrdcp -f cmsswEvents.root root://eoscms.cern.ch/${eosOutputPath}/cmsswEvents_Run${run}.root
