import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('runNumber',
                496,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'RunNumber.'
                )

options.register('beamEnergy',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Beam energy.'
                )

options.register('beamParticlePDGID',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Beam particles PDG ID.'
                )

options.register('runType',
                 "Beam",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Run type: Pedestal, Beam, Simulation.'
                )

options.register('setupConfiguration',
                5,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'setupConfiguration (1: July - 4: 20 Layers in October in H6A".'
                )

options.register('NHexaBoards',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('SubtractPedestal',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Subtract the pedestals.'
                )

options.register('electronicMap',
                 'map_DESY_March2018_config4_V0.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('adcCalibrations',
                 'hgcal_calibration.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal ADC to MIP calibration file in HGCal/CondObjects/data/')

options.register('hgcalLayout',
                 'layerGeom_desymarch2018_configuration4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')


options.register('layerPositionFile',
                 '/afs/cern.ch/user/t/tquast/CMSSW_8_3_0/src/HGCal/CondObjects/data/layer_distances_DESY_March2018_config4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 ''
                )

options.maxEvents = -1

options.parseArguments()
print options


electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap
hgcalLayout="HGCal/CondObjects/data/%s" % options.hgcalLayout
adcCalibrations="HGCal/CondObjects/data/%s" % options.adcCalibrations

################################
process = cms.Process("analysis")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################

#only read time stamps for beam runs
readTimeStamps = False
if options.beamParticlePDGID==0:
    readTimeStamps = False




process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/ORM_raw/HexaData_Run%04d.raw"%options.runNumber),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30787),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(12),          
                            NumberOfBytesForTheTrailer=cms.untracked.uint32(4),         
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(12),  
                            NSkipEvents=cms.untracked.uint32(0),
                            ReadTimeStamps=cms.untracked.bool(readTimeStamps),
                            DataFormats=cms.untracked.uint32(1),
                            timingFiles=cms.vstring(),
                            beamEnergy=cms.untracked.double(options.beamEnergy),
                            beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),
                            runType=cms.untracked.string(options.runType),
                            setupConfiguration=cms.untracked.uint32(options.setupConfiguration)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("/home/tquast/tbMarch2018_DESY/analysis/analysed_%04d.root"%options.runNumber))


process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        GlobalTimestampCollectionName=cms.string("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicMap=cms.untracked.string(electronicMap),
                                        SubtractPedestal=cms.untracked.bool(bool(1)),
                                        MaskNoisyChannels=cms.untracked.bool(bool(0)),
                                        HighGainPedestalFileName=cms.untracked.string("/home/tquast/tbMarch2018_DESY/pedestals/pedestalHG_%04d.txt"%(options.runNumber)),
                                        LowGainPedestalFileName=cms.untracked.string("/home/tquast/tbMarch2018_DESY//pedestals/pedestalLG_%04d.txt"%(options.runNumber)),
                                        NoisyChannelsFileName=cms.untracked.string("/home/tquast/tbMarch2018_DESY/pedestals/noisyChannels_%04d.txt"%(options.runNumber))
)

process.rawhitplotter = cms.EDAnalyzer("RawHitPlotter",
                                       InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       SubtractCommonMode=cms.untracked.bool(True)
)

process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                        RUNDATA = cms.InputTag("source", "RunData" ), 
                                        GlobalTimestampCollectionName=cms.InputTag("rawhitproducer","HGCALGLOBALTIMESTAMPS"),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        MaskNoisyChannels=cms.untracked.bool(bool(0)),
                                        ChannelsToMaskFileName=cms.untracked.string("/home/tquast/pedestals/noisyChannels_%04d.txt"%(options.runNumber)),
                                        NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                        TimeSample3ADCCut = cms.untracked.double(15.),
                                        investigatePulseShape = cms.untracked.bool(True),
                                        timingNetworks = cms.untracked.string("")
)

process.rechitplotter = cms.EDAnalyzer("RecHitPlotter",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       MipThreshold=cms.untracked.double(200),
                                       NoiseThreshold=cms.untracked.double(20)
)


process.cellenergyplotting = cms.EDAnalyzer("CellEnergyPlotter",
                                RUNDATA = cms.InputTag("source", "RunData" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                HGCALTBCOMMONMODENOISE = cms.InputTag("rechitproducer","HGCALTBCOMMONMODENOISEMAP" ),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                commonModeNoiseRejectionType = cms.int32(0)       #0: none, else 1-..., default: 0 
)

process.variablecomputation = cms.EDProducer("VariableComputation",
                                RUNDATA = cms.InputTag("source", "RunData" ), 
                                MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                DWCTRACKS = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                UserRecordCollectionName=cms.untracked.string("VariableUserRecords"),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                layerPositionFile=cms.string(options.layerPositionFile),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                NLayers=cms.untracked.int32(options.NHexaBoards),
                                DNNInputFile = cms.untracked.string(""),
                                NColorsInputImage = cms.untracked.int32(-1),
                                CellEnergyCut = cms.untracked.double(0.5)
)

process.ntupelizer = cms.EDAnalyzer("NTupelizer",
                                USERRECORDS = cms.InputTag("variablecomputation","VariableUserRecords" ),
                                UserRecordKeys = cms.vstring(["NRechits", "I_EV1", "I_EV2", "E1_tot", "E7_tot", "E19_tot", "E19_tot", "E37_tot", "EAll_tot", "depthX0", "EAll_layer1", "EAll_layer2", "EAll_layer3", "NAll_layer1", "NAll_layer2", "NAll_layer3", "E1PerE7_layer1", "E1PerE7_layer2", "E1PerE7_layer3", "E7PerE19_layer1", "E7PerE19_layer2", "E7PerE19_layer3"])
)



####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


#process.p = cms.Path(process.rawhitproducer*process.rawhitplotter*process.rechitproducer*process.rechitplotter*process.cellenergyplotting*process.variablecomputation*process.ntupelizer)
process.p = cms.Path(process.rawhitproducer*process.rawhitplotter*process.rechitproducer*process.rechitplotter*process.cellenergyplotting*process.variablecomputation*process.ntupelizer)

