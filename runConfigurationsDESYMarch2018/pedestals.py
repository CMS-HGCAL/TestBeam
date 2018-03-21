import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
options.register('dataFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/ORM_raw/HexaData_Run0492.raw',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'file containing raw input')

options.register('outputFile',
                 '/afs/cern.ch/user/t/tquast/Desktop/tb2018_DESY/quickAnalysis/pedestals/pedestal_0492.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('pedestalHighGainFile',
                 '/afs/cern.ch/user/t/tquast/Desktop/tb2018_DESY/quickAnalysis/pedestals/pedestalHG_0492.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('pedestalLowGainFile',
                 '/afs/cern.ch/user/t/tquast/Desktop/tb2018_DESY/quickAnalysis/pedestals/pedestalLG_0492.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('noisyChannelsFile',
                 '/afs/cern.ch/user/t/tquast/Desktop/tb2018_DESY/quickAnalysis/pedestals/noisyChannels_0492.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('electronicMap',
                 'map_CERN_Hexaboard_July_6Layers.txt',
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
                4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )



options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
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
                            fileNames=cms.untracked.vstring("file:%s"%(options.dataFile)),
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

process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))
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
                                         HighGainPedestalFileName=cms.untracked.string(options.pedestalHighGainFile),
                                         LowGainPedestalFileName=cms.untracked.string(options.pedestalLowGainFile),
                                         WriteNoisyChannelsFile=cms.untracked.bool(True),
                                         NoisyChannelsFileName=cms.untracked.string(options.noisyChannelsFile),
)


process.p = cms.Path( process.rawdataplotter*process.pedestalplotter )


