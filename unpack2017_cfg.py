import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/disk2_2TB/July2017_TB_data_orm',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('fileName',
                 'Module63.raw',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Input file to process')

options.register('outputFolder',
                 '/afs/cern.ch/work/a/asteen/public/data/july2017/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('compressedData',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Option to set if the data have beem compressed')

options.register('electronicMap',
                 "HGCal/CondObjects/data/map_CERN_Hexaboard_September_17Sensors_7EELayers_10FHLayers_V1.txt",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')


options.maxEvents = -1
options.output = "cmsswEvents.root"

options.parseArguments()
print options
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

pedestalHighGain=options.outputFolder+"/pedestalHG.txt"
pedestalLowGain=options.outputFolder+"/pedestalLG.txt"
noisyChannels=options.outputFolder+"/noisyChannels.txt"

################################
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################

numberOfBytesForTheHeader=48
numberOfBytesForTheEventTrailers=2
numberOfBytesPerReadOut=30784
if options.compressedData==1:
    numberOfBytesPerReadOut=15392
process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string(options.electronicMap),
                            fileNames=cms.untracked.vstring("file:%s/%s"%(options.dataFolder,options.fileName)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOfBytesPerReadOut=cms.untracked.uint32(numberOfBytesPerReadOut),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(numberOfBytesForTheHeader),
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(numberOfBytesForTheEventTrailers),
                            NSkipEvents=cms.untracked.uint32(0),
                            CompressedData=cms.untracked.bool(options.compressedData)
)

filename = options.outputFolder+"/PedestalOutput.root"
process.TFileService = cms.Service("TFileService", fileName=cms.string(filename))

process.output = cms.OutputModule("PoolOutputModule",fileName = cms.untracked.string(options.output))

process.content = cms.EDAnalyzer("EventContentAnalyzer") #add process.content in cms.Path if you want to check which collections are in the event

process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(True),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(options.electronicMap),
                                         HighGainPedestalFileName=cms.untracked.string(pedestalHighGain),
                                         LowGainPedestalFileName=cms.untracked.string(pedestalLowGain),
                                         WriteNoisyChannelsFile=cms.untracked.bool(True),
                                         NoisyChannelsFileName=cms.untracked.string(noisyChannels),
                                         NTSForPedestalComputation=cms.untracked.int32(2),

)

process.rawdataplotter = cms.EDAnalyzer("RawDataPlotter",
                                        SensorSize=cms.untracked.int32(128),
                                        EventPlotter=cms.untracked.bool(False),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata")
)

process.p = cms.Path( process.pedestalplotter )

process.end = cms.EndPath(process.output)

