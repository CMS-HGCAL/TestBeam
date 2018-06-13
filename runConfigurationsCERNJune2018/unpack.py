import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:


options.register('dataFile',
                 '/home/tquast/tbJune2018_H2/EUDAQ_unpacked/run000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 223,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('timingFile',
                 '/home/tquast/tbJune2018_H2/timingFiles/timing000223.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('outputFile',
                 '/home/tquast/tbJune2018_H2/pedestals/pedestals_000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('unpackedFile',
                 '/home/tquast/tbJune2018_H2/CMSSW_unpacked/unpacked_000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('beamEnergy',
                120,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Beam energy.'
                )

options.register('beamParticlePDGID',
                13,
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
                18,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'setupConfiguration (1: July - 4: 20 Layers in October in H6A".'
                )

options.register('pedestalHighGainFile',
                 '/home/tquast/tbJune2018_H2/pedestals/pedestalHG_000223.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('pedestalLowGainFile',
                 '/home/tquast/tbJune2018_H2/pedestals/pedestalLG_000223.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('noisyChannelsFile',
                 '/home/tquast/tbJune2018_H2/pedestals/noisyChannels_000223.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('electronicMap',
                 "map_CERN_Hexaboard_June_28Sensors_28EELayers_V0.txt",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')

options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '.'
                )


options.maxEvents = -1

options.parseArguments()
print options

electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap

################################
process = cms.Process("CMSSWUnpacker")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

################################
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################

process.source = cms.Source("HGCalTBEUDAQDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:%s"%options.dataFile),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NSkipEvents=cms.untracked.uint32(0),
                            runNumber=cms.untracked.int32(options.runNumber),
                            beamEnergy=cms.untracked.double(options.beamEnergy),
                            beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),
                            runType=cms.untracked.string(options.runType),
                            setupConfiguration=cms.untracked.uint32(options.setupConfiguration)
)

process.timingfilewriter = cms.EDAnalyzer("HGCalTBTimingFileWriter",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        TimingFilePath=cms.untracked.string(options.timingFile)
)

process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(True),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(electronicMap),
                                         HighGainPedestalFileName=cms.untracked.string(options.pedestalHighGainFile),
                                         LowGainPedestalFileName=cms.untracked.string(options.pedestalLowGainFile),
                                         WriteNoisyChannelsFile=cms.untracked.bool(True),
                                         NoisyChannelsFileName=cms.untracked.string(options.noisyChannelsFile),
)

process.p = cms.Path( process.timingfilewriter * process.pedestalplotter)

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.unpackedFile)
)

process.end = cms.EndPath(process.output)

