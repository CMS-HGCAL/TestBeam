import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/eos/user/t/tquast/data/Testbeam/May2017',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input.')

options.register('runNumber',
                 151,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process.')

options.register('beamMomentum',
                 250.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Momentum of the electron beam.')

options.register('outputPath',
                 '/eos/user/t/tquast/outputs/Testbeam/May2017/LayerSumAnalyzer/run_xx.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored.')

options.register('reportEvery',
                100,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                'Frequency of event count printouts on the console.')


options.maxEvents = -1

options.parseArguments()
print options
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

################################
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.load('HGCal.StandardSequences.LocalReco_cff')

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

####################################

process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt"),
                            fileNames=cms.untracked.vstring("file:%s/HexaData_Run%04d.raw"%(options.dataFolder,options.runNumber)),#the current version can handle only one file
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NOrmBoards=cms.untracked.uint32(1),
                            NHexaBoards=cms.untracked.uint32(1),
                            NSkirocsPerHexa=cms.untracked.uint32(4),
                            NChannelsPerSkiroc=cms.untracked.uint32(64),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30788)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputPath))


process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata")
                                        )


process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag('rawhitproducer','HGCALTBRAWHITS'),
                                        LowGainPedestalFileName = cms.string("/afs/cern.ch/user/t/tquast/CMSSW_8_0_0_pre5/src/HGCal/CondObjects/data/pedestal_LG_OneLayer2017.txt"),
                                        HighGainPedestalFileName = cms.string('/afs/cern.ch/user/t/tquast/CMSSW_8_0_0_pre5/src/HGCal/CondObjects/data/pedestal_HG_OneLayer2017.txt'),
                                        LG2HG = cms.untracked.vdouble(10.0),
                                        TOT2LG = cms.untracked.vdouble(10.0),
                                        HighGainADCSaturation = cms.untracked.double(1800),
                                        LowGainADCSaturation = cms.untracked.double(1800),
                                        ElectronicsMap = cms.untracked.string('HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt'),
                                        CommonModeThreshold = cms.untracked.double(3.),
                                        KeepOnlyTimeSample3 = cms.untracked.bool(True),
                                        performParabolicFit = cms.untracked.bool(False)
                                        )


process.LayerSumAnalyzer2017 = cms.EDAnalyzer("Layer_Sum_Analyzer_2017",
                                  HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                  layers_config = cms.int32(1),
                                  mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt'),
                                  run = cms.int32(options.runNumber),
                                  beamMomentum = cms.double(options.beamMomentum)
                              )


process.p = cms.Path( process.rawhitproducer * process.rechitproducer * process.LayerSumAnalyzer2017)

#process.p = cms.Path( process.rawhitproducer * process.hgcaltbrechits2017 * process.LayerSumAnalyzer2017)
