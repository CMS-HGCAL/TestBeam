import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/afs/cern.ch/work/b/barneyd/public/data_May_2017/disk2_2TB/eudaq_data',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 850,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('outputFolder',
                 '/afs/cern.ch/work/a/asteen/public/data/may2017/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.maxEvents = -1
options.output = "cmsswEvents.root"

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

process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/HexaOutput_%d.root"%(options.outputFolder,options.runNumber)))

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output)
                                  )


pedestalHighGain="pedestalHG_125.txt"
pedestalLowGain="pedestalLG_125.txt"

process.rawdataplotter = cms.EDAnalyzer("RawDataPlotter",
                                        SensorSize=cms.untracked.int32(128),
                                        EventPlotter=cms.untracked.bool(False),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
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
                                       SubtractPedestal=cms.untracked.bool(False),
                                       SubtractCommonMode=cms.untracked.bool(False),
                                       CommonModeThreshold=cms.untracked.double(3),
                                       HighGainPedestalFileName=cms.string(pedestalHighGain),
                                       LowGainPedestalFileName=cms.string(pedestalLowGain)
                                       )

process.pulseshapeplotter = cms.EDAnalyzer("PulseShapePlotter",
                                           InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                           ElectronicMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt"),
                                           HighGainPedestalFileName=cms.string(pedestalHighGain),
                                           LowGainPedestalFileName=cms.string(pedestalLowGain)
                                           )


process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag('rawhitproducer','HGCALTBRAWHITS'),
                                        LowGainPedestalFileName = cms.string('pedestalLG_125.txt'),
                                        HighGainPedestalFileName = cms.string('pedestalHG_125.txt'),
                                        LG2HG = cms.untracked.vdouble(10.0),
                                        TOT2LG = cms.untracked.vdouble(10.0),
                                        HighGainADCSaturation = cms.untracked.double(1800),
                                        LowGainADCSaturation = cms.untracked.double(1800),
                                        ElectronicsMap = cms.untracked.string('HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt'),
                                        CommonModeThreshold = cms.untracked.double(3.),
                                        KeepOnlyTimeSample3 = cms.untracked.bool(True)
                                        )

process.rechitplotter = cms.EDAnalyzer("RecHitPlotter",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       ElectronicMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt"),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       NoiseThreshold=cms.untracked.double(15)
                                       )

#process.p = cms.Path( process.rawdataplotter )
#process.p = cms.Path( process.rawhitproducer*process.rawhitplotter )
#process.p = cms.Path( process.rawhitproducer*process.rawhitplotter*process.pulseshapeplotter )
process.p = cms.Path( process.rawhitproducer*process.rechitproducer*process.rechitplotter )

process.end = cms.EndPath(process.output)

