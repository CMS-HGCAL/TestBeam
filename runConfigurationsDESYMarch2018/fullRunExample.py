import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

#full run configuration for DESY (March 2018) run 1174
#author: Thorben Quast, thorben.quast@cern.ch
#18 May 2018

#execution on any machine with lxplus 
#1. source /cvmfs/cms.cern.ch/cmsset_default.sh;
#2. cd <CMSSW-dir>/src/<Testbeam-dir>
#3. cmsenv
#4. cmsRun runConfigurationsDESYMarch2018/fullRunExample.py

############################hard coded configuration: adjust for each run##########################
runNumber = 1174
beamEnergy = 6.0
setupConfiguration = 6
SkipFirstNEventsInTelescopeFile = 1

############################hard coded configuration: fixed for all runs documented in the according Google spreadsheet from 28th March 2018#########################
ORM_dataFolder = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/ORM_raw"
telescope_dataFolder = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/beamTelescopeReco_corryvreckan/v1"
beamParticlePDGID = 11
runType = "Beam"

pedestalHighGainFile = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/pedestals/v1/pedestalsHG_%s.txt" % runNumber
pedestalLowGainFile = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/pedestals/v1/pedestalsLG_%s.txt" % runNumber
noisyChannelsFile = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/pedestals/v1/noisyChannels_%s.txt" % runNumber
VariablesToPlot = ["xmean", "ymean", "NRechits", "E7_tot", "E19_tot", "EAll_tot", "EAllHG_tot", "EAllLG_tot", "EAllTOT_tot", "EAll_layer1", "EAll_layer2", "EAll_layer3" ]
eventNumbersToDisplay = range(1, 11)        #displays for the first ten events
outputFile = "output_%s.root" % runNumber


NHexaBoards = 3
electronicMap = "HGCal/CondObjects/data/map_CERN_Hexaboard_DESY_3Sensors_3EELayers_V1.txt"
hgcalLayout = "HGCal/CondObjects/data/layerGeom_desymarch2018_configuration4.txt"
adcCalibrations = "HGCal/CondObjects/data/hgcal_calibration_DESYMarch2018.txt"
layerPositions = "/afs/cern.ch/user/t/tquast/CMSSW_9_3_0/src/HGCal/CondObjects/data/layer_distances_DESY_March2018_config4.txt"     #needs to be the absolute path
PIStagePositionFile = "/afs/cern.ch/user/t/tquast/CMSSW_9_3_0/src/HGCal/CondObjects/data/pi_stage_positions_DESYMarch2018.txt"     #needs to be the absolute path

reportEvery = 1000
maxEvents = -1 #-1 for analysing all events                 
#############################################################################



################################
process = cms.Process("analysis")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery
####################################


process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:%s/HexaData_Run%04d.raw"%(ORM_dataFolder, runNumber)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30787),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(12),
                            NumberOfBytesForTheTrailer=cms.untracked.uint32(4),
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(12),  
                            NSkipEvents=cms.untracked.uint32(0),
                            ReadTimeStamps=cms.untracked.bool(True),
                            DataFormats=cms.untracked.uint32(1),
                            timingFiles=cms.vstring(),
                            beamEnergy=cms.untracked.double(beamEnergy),
                            beamParticlePDGID=cms.untracked.int32(beamParticlePDGID),
                            runType=cms.untracked.string(runType),
                            setupConfiguration=cms.untracked.uint32(setupConfiguration)
)

process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        GlobalTimestampCollectionName=cms.string("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicMap=cms.untracked.string(electronicMap),
                                        SubtractPedestal=cms.untracked.bool(True),
                                        MaskNoisyChannels=cms.untracked.bool(True),
                                        HighGainPedestalFileName=cms.untracked.string(pedestalHighGainFile),
                                        LowGainPedestalFileName=cms.untracked.string(pedestalLowGainFile),
                                        ChannelsToMaskFileName=cms.untracked.string(noisyChannelsFile)
)

process.rawhitplotter = cms.EDAnalyzer("RawHitPlotter",
                                       InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       NHexaBoards=cms.untracked.int32(NHexaBoards),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       SubtractCommonMode=cms.untracked.bool(True)
)

process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag("rawhitproducer", "HGCALTBRAWHITS"),
                                        RUNDATA = cms.InputTag("source", "RunData" ), 
                                        GlobalTimestampCollectionName=cms.InputTag("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        MaskNoisyChannels=cms.untracked.bool(bool(0)),
                                        ChannelsToMaskFileName=cms.untracked.string(""),
                                        NHexaBoards=cms.untracked.int32(NHexaBoards),
                                        TimeSample3ADCCut = cms.untracked.double(15.),
                                        investigatePulseShape = cms.untracked.bool(True),
                                        timingNetworks = cms.untracked.string("")
)

process.rechitplotter = cms.EDAnalyzer("RecHitPlotter",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       NHexaBoards=cms.untracked.int32(NHexaBoards),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       MipThreshold=cms.untracked.double(200),
                                       NoiseThreshold=cms.untracked.double(20)
)


process.daturaproducer = cms.EDProducer("HGCalTBDATURATelescopeProducer",
                            OutputCollectionName = cms.string("HGCalTBDATURATracks"), 
                            RUNDATA = cms.InputTag("source","RunData"),
                            inputFile = cms.string("%s/reco_%06d.root"%(telescope_dataFolder, runNumber)),
                            SkipFirstNEventsInTelescopeFile = cms.int32(SkipFirstNEventsInTelescopeFile),  #first event in the beam telescope file is a BORE event
                            layerPositionFile = cms.string(layerPositions),
                            PIStagePositionFile = cms.string(PIStagePositionFile)
)

process.cellenergyplotting = cms.EDAnalyzer("CellEnergyPlotter",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                HGCALTBCOMMONMODENOISE = cms.InputTag("rechitproducer","HGCALTBCOMMONMODENOISEMAP" ),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                NHexaBoards=cms.untracked.int32(NHexaBoards),
                                commonModeNoiseRejectionType = cms.int32(0)       #0: none, else 1-..., default: 0 
)

process.correlationanalysis = cms.EDAnalyzer("DATURATelescopeCorrelator",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                DATURATelescopeData = cms.InputTag("daturaproducer","HGCalTBDATURATracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                NHexaBoards=cms.untracked.int32(NHexaBoards),
                                n_bins = cms.int32(50),
                                max_dim_x_DUT = cms.double(10.),
                                max_dim_y_DUT = cms.double(10.)
)

process.variablecomputation = cms.EDProducer("VariableComputation",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                DWCTRACKS = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                UserRecordCollectionName=cms.untracked.string("VariableUserRecords"),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                layerPositionFile=cms.string(layerPositions),
                                NHexaBoards=cms.untracked.int32(NHexaBoards),
                                NLayers=cms.untracked.int32(NHexaBoards),
                                NColorsInputImage = cms.untracked.int32(-1),
                                CellEnergyCut = cms.untracked.double(0.5)
)

process.variablewriter = cms.EDAnalyzer("NTupelizer",
                                USERRECORDS = cms.InputTag("variablecomputation","VariableUserRecords" ),
                                UserRecordKeys = cms.vstring(VariablesToPlot)
)

process.eventdisplay = cms.EDAnalyzer("EventDisplay",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                MWCHAMBERS = cms.InputTag("","" ), 
                                DWCTRACKS = cms.InputTag("","" ), 
                                electronicsMap = cms.untracked.string(electronicMap),
                                NHexaBoards=cms.untracked.int32(NHexaBoards),
                                eventsToPlot=cms.vint32(eventNumbersToDisplay)
                              )

process.rechitntupler = cms.EDAnalyzer("RecHitNtupler",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       layerPositionFile = cms.untracked.string(layerPositions),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(True),
                                       MipThreshold=cms.untracked.double(5.0),
                                       NoiseThreshold=cms.untracked.double(0.5)
)

process.trackimpactntupler = cms.EDAnalyzer("ImpactPointNtupler",
                                       DATURATelescopeData = cms.InputTag("daturaproducer","HGCalTBDATURATracks" ),
                                       RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ),
                                       nLayers=cms.untracked.int32(NHexaBoards)
)

####################################


process.p = cms.Path(process.rawhitproducer*process.rawhitplotter*process.rechitproducer*process.rechitplotter*process.daturaproducer*process.cellenergyplotting*process.correlationanalysis*process.variablecomputation*process.variablewriter*process.eventdisplay*process.rechitntupler*process.trackimpactntupler)


