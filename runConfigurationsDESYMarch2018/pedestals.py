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


options.register('electronicMap',
                 'map_DESY_March2018_config4_V0.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('NTSForPedestalComputation',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of time samples used for pedestal computation.'
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

options.register('dataFormat',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Data formats int -> important for parameter setting')


options.register('NHexaBoards',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )


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


process = cms.Process("pedestals")

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

####################################

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


#only read time stamps for beam runs
readTimeStamps = False
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
                            fileNames=cms.untracked.vstring("file:/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/ORM_raw/HexaData_Run%04d.raw"%(options.runNumber)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30787),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(numberOfBytesForTheHeader),          #for the new headers/trailers from run 1241 onward: 12
                            NumberOfBytesForTheTrailer=cms.untracked.uint32(4),         #for the new headers/trailers from run 1241 onward: 4
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(numberOfBytesForTheEventTrailers),   #for the new headers/trailers from run 1241 onward: 12
                            NSkipEvents=cms.untracked.uint32(0),
                            ReadTimeStamps=cms.untracked.bool(readTimeStamps),
                            DataFormats=cms.untracked.uint32(readTimeStampFormat),
                            timingFiles=cms.vstring(),
                            beamEnergy=cms.untracked.double(options.beamEnergy),
                            beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),
                            runType=cms.untracked.string(options.runType),
                            setupConfiguration=cms.untracked.uint32(options.setupConfiguration)
)


################################
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################


process.content = cms.EDAnalyzer("EventContentAnalyzer") #add process.content in cms.Path if you want to check which collections are in the event

process.TFileService = cms.Service("TFileService", fileName=cms.string('/home/tquast/pedestals/pedestal_%04d.root'%(options.runNumber)))
process.rawdataplotter = cms.EDAnalyzer("RawDataPlotter",
                                        SensorSize=cms.untracked.int32(128),
                                        EventPlotter=cms.untracked.bool(False),
                                        NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata")
)

process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(True),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(electronicMap),
                                         NTSForPedestalComputation=cms.untracked.int32(options.NTSForPedestalComputation),
                                         HighGainPedestalFileName=cms.untracked.string("/home/tquast/tbMarch2018_DESY/pedestals/pedestalHG_%04d.txt"%(options.runNumber)),
                                         LowGainPedestalFileName=cms.untracked.string("/home/tquast/tbMarch2018_DESY/pedestals/pedestalLG_%04d.txt"%(options.runNumber)),
                                         WriteNoisyChannelsFile=cms.untracked.bool(True),
                                         NoisyChannelsFileName=cms.untracked.string("/home/tquast/tbMarch2018_DESY/pedestals/noisyChannels_%04d.txt"%(options.runNumber)),
)


process.p = cms.Path( process.rawdataplotter*process.pedestalplotter )


