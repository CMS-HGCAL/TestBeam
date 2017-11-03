import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/july2017/July2017_TB_data_orm/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 106,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('outputFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/skiroc2cmsdata/skiroc2cms.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where unpacked data is stored')

options.register('electronicMap',
                 'map_CERN_Hexaboard_July_6Layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('beamEnergy',
                250,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Beam energy.'
                )

options.register('beamParticlePDGID',
                211,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Beam particles PDG ID.'
                )

options.register('setupConfiguration',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'setupConfiguration (1: July - 4: 20 Layers in October in H6A".'
                )

options.register('dataFormat',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Data formats int -> important for parameter setting')

options.register('reportEvery',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '.'
                )

options.maxEvents = -1

options.parseArguments()
print options
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap

################################
process = cms.Process("Unpack")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

####################################

#only read time stamps for beam runs
readTimeStamps = True
if options.beamParticlePDGID==0:
    readTimeStamps = False

if options.dataFormat==0 :
    numberOfBytesForTheHeader=8
    numberOfBytesForTheTrailer=4
    numberOfBytesForTheEventTrailers=4
    readTimeStampFormat = 0     #from the txt file
elif options.dataFormat==1 :
    numberOfBytesForTheHeader=12
    numberOfBytesForTheTrailer=4
    numberOfBytesForTheEventTrailers=12
    readTimeStampFormat = 1 #from the raw file


process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:%s/HexaData_Run%04d.raw"%(options.dataFolder,options.runNumber)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30787),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(numberOfBytesForTheHeader),          #for the new headers/trailers from run 1241 onward: 12
                            NumberOfBytesForTheTrailer=cms.untracked.uint32(4),         #for the new headers/trailers from run 1241 onward: 4
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(numberOfBytesForTheEventTrailers),   #for the new headers/trailers from run 1241 onward: 12
                            NSkipEvents=cms.untracked.uint32(0),
                            ReadTimeStamps=cms.untracked.bool(readTimeStamps),
                            DataFormats=cms.untracked.uint32(readTimeStampFormat),
                            timingFiles=cms.vstring("%s/HexaData_Run%04d_TIMING_RDOUT_ORM0.txt"%(options.dataFolder,options.runNumber),
                                                    "%s/HexaData_Run%04d_TIMING_RDOUT_ORM1.txt"%(options.dataFolder,options.runNumber),
                                                    "%s/HexaData_Run%04d_TIMING_RDOUT_ORM2.txt"%(options.dataFolder,options.runNumber)),
                            beamEnergy=cms.untracked.uint32(options.beamEnergy),
                            beamParticlePDGID=cms.untracked.double(options.beamParticlePDGID),
                            setupConfiguration=cms.untracked.uint32(options.setupConfiguration)
)


process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.outputFile)
)

process.p = cms.Path( )

process.end = cms.EndPath(process.output)

