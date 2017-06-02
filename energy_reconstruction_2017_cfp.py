import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/eos/user/t/tquast/data/Testbeam/May2017',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 151,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('beamMomentum',
                 250.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Momentum of the electron beam.')

options.register('outputFolder',
                 '/eos/user/t/tquast/outputs/Testbeam/May2017',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

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

process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/EnergySums_%04d.root"%(options.outputFolder,options.runNumber)))



process.content = cms.EDAnalyzer("EventContentAnalyzer") #add process.content in cms.Path if you want to check which collections are in the event

process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata")
                                        )


process.LayerSumAnalyzer2017 = cms.EDAnalyzer("Layer_Sum_Analyzer_2017",
                                  HGCALTBRECHITS = cms.InputTag("hgcaltbrechits2017","","unpack" ),
                                  layers_config = cms.int32(1),
                                  mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt'),
                                  run = cms.int32(options.runNumber),
                                  beamMomentum = cms.double(options.beamMomentum)
                              )


process.p = cms.Path( process.rawhitproducer * process.hgcaltbrechits2017 * process.LayerSumAnalyzer2017)

