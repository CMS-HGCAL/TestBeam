import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/afs/cern.ch/work/a/asteen/public/data/may2017/raw',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 850,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input file to process')

options.register('outputFolder',
                 '/tmp/asteen',
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
process.source = cms.Source("HGCalTBRawDataSource",
                            #FilePath=cms.untracked.string(options.dataFolder),
                            #FileName=cms.untracked.string("file:HexaData_Run00%d.raw"%(options.runNumber)),
                            ElectronicMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt"),
                            fileNames=cms.untracked.vstring("file:%s/HexaData_Run00%d.raw"%(options.dataFolder,options.runNumber)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NOrmBoards=cms.untracked.uint32(1),
                            NHexaBoards=cms.untracked.uint32(1),
                            NSkirocsPerHexa=cms.untracked.uint32(4),
                            NChannelsPerSkiroc=cms.untracked.uint32(64),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30788)
)

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output)
                                  )

process.end = cms.EndPath(process.output)

#help(process)
