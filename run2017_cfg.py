import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 #'./',
								 #'/afs/cern.ch/work/a/asteen/public/data/may2017/raw',
                 '/afs/cern.ch/work/b/barneyd/public/data_May_2017/disk2_2TB/eudaq_data',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 850,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input file to process')

options.register('outputFolder',
                 '/afs/cern.ch/work/a/asteen/public/data/may2017/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Result of processing')

options.maxEvents = -1
options.output = "toto.root"

options.parseArguments()
print options
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

################################
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################
#process.load('HGCal.StandardSequences.RawToDigi_cff')

process.source = cms.Source("HGCalTBRawDataSource",
                            #FilePath=cms.untracked.string(options.dataFolder),
                            #FileName=cms.untracked.string("file:HexaData_Run00%d.raw"%(options.runNumber)),
                            ElectronicMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt"),
                            fileNames=cms.untracked.vstring("file:%s/HexaData_Run%04d.raw"%(options.dataFolder,options.runNumber)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NOrmBoards=cms.untracked.uint32(1),
                            NHexaBoards=cms.untracked.uint32(1),
                            NSkirocsPerHexa=cms.untracked.uint32(4),
                            NChannelsPerSkiroc=cms.untracked.uint32(64),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30788)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/toto_%d.root"%(options.outputFolder,options.runNumber)))

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output)
                                  )


pedestalHighGain="pedestalHG_125.txt"
pedestalLowGain="pedestalLG_125.txt"

process.rawdataplotter = cms.EDAnalyzer("RawDataPlotter",
                                        SensorSize=cms.untracked.int32(128),
                                        EventPlotter=cms.untracked.bool(False),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        HighGainPedestalFileName=cms.string(pedestalHighGain),
                                        LowGainPedestalFileName=cms.string(pedestalLowGain)
                                        )
process.content = cms.EDAnalyzer("EventContentAnalyzer") #add process.content in cms.Path if you want to check which collections are in the event

process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata")
                                        )

process.rawhitplotter = cms.EDAnalyzer("RawHitPlotter",
                                       InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                       ElectronicMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt"),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       SubtractPedestal=cms.untracked.bool(True),
                                       SubtractCommonMode=cms.untracked.bool(True),
                                       HighGainPedestalFileName=cms.string(pedestalHighGain),
                                       LowGainPedestalFileName=cms.string(pedestalLowGain)
                                       )




#process.p = cms.Path( process.rawdataplotter )
process.p = cms.Path( process.rawhitproducer*process.rawhitplotter )

process.end = cms.EndPath(process.output)

#help(process)
