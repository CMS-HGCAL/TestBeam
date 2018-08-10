import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:

options.register('fileName',
                 './rawdata.raw',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Input file (path needed) to process')

options.register('outputFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('compressedData',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Option to set if the data have beem compressed')

options.register('electronicMap',
                 "HGCal/CondObjects/data/map_CERN_Hexaboard_OneModule_V2.txt",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')

options.register('writePedestal',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'set to true if pedestal files have to be created')

options.register('writeNoise',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'set to true if noise files have to be created')

options.maxEvents = -1

options.parseArguments()

aname=os.path.basename(options.fileName)
aname,extent=os.path.splitext(aname)
aname=aname+"_edm.root"
options.output=options.outputFolder+aname

print options

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
                            fileNames=cms.untracked.vstring("file:%s"%(options.fileName)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOfBytesPerReadOut=cms.untracked.uint32(numberOfBytesPerReadOut),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(numberOfBytesForTheHeader),
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(numberOfBytesForTheEventTrailers),
                            NSkipEvents=cms.untracked.uint32(0),
                            CompressedData=cms.untracked.bool(options.compressedData)
)

outputFileName=os.path.basename(options.fileName)
outputFileName,extent=os.path.splitext(outputFileName)
outputFileName=options.outputFolder+outputFileName+"_pedestal.root"
process.TFileService = cms.Service("TFileService", fileName=cms.string(outputFileName))

process.output = cms.OutputModule("PoolOutputModule",fileName = cms.untracked.string(options.output))

process.content = cms.EDAnalyzer("EventContentAnalyzer") #add process.content in cms.Path if you want to check which collections are in the event


process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(options.writePedestal),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(options.electronicMap),
                                         HighGainPedestalFileName=cms.untracked.string(pedestalHighGain),
                                         LowGainPedestalFileName=cms.untracked.string(pedestalLowGain),
                                         WriteNoisyChannelsFile=cms.untracked.bool(options.writeNoise),
                                         NoisyChannelsFileName=cms.untracked.string(noisyChannels),
                                         NTSForPedestalComputation=cms.untracked.int32(7),

)

process.treeproducer = cms.EDAnalyzer("TreeProducer",
                                      InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                      SubtractPedestal=cms.untracked.bool(False),
                                      HighGainPedestalFileName=cms.untracked.string(pedestalHighGain),
                                      LowGainPedestalFileName=cms.untracked.string(pedestalLowGain),
)

process.p = cms.Path( process.pedestalplotter*process.treeproducer )

process.end = cms.EndPath(process.output)

