import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:


options.register('dataFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/skiroc2cmsdata/SKIROC2_CMS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('processedFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/rawhits/RAWHITS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('electronicMap',
                 'map_CERN_Hexaboard_July_6Layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('NHexaBoards',
                10,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('pedestalHighGainFile',
                 '/home/tquast/tb2017/pedestals/pedestalHG.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('pedestalLowGainFile',
                 '/home/tquast/tb2017/pedestals/pedestalLG.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('noisyChannelsFile',
                 '/home/tquast/tb2017/pedestals/noisyChannels.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('SubtractPedestal',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Subtract the pedestals.'
                )

options.register('MaskNoisyChannels',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.register('reportEvery',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.maxEvents = -1

options.parseArguments()
print options


electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap


################################
process = cms.Process("rawhitproducer")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s"%options.dataFile)
)


process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.processedFile),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRAWHITS_*',
                                                                         'keep *_*_HGCALGLOBALTIMESTAMPS_*',
                                                                         'keep *_*_RunData_*')
)



process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        GlobalTimestampCollectionName=cms.string("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicMap=cms.untracked.string(electronicMap),
                                        SubtractPedestal=cms.untracked.bool(bool(options.SubtractPedestal)),
                                        MaskNoisyChannels=cms.untracked.bool(bool(options.MaskNoisyChannels)),
                                        HighGainPedestalFileName=cms.untracked.string(options.pedestalHighGainFile),
                                        LowGainPedestalFileName=cms.untracked.string(options.pedestalLowGainFile),
                                        ChannelsToMaskFileName=cms.untracked.string(options.noisyChannelsFile)
)


process.p = cms.Path( process.rawhitproducer )

process.end = cms.EndPath(process.output)
