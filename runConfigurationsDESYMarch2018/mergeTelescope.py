import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:


options.register('dataFile',
                 '/home/tquast/tbMarch2018_DESY/reco/RECHITS_1174.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('TelescopeFile',
                 '/home/tquast/tbMarch2018_DESY/beamTelescope/reco_001174.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the reconstructed DATURA file.')

options.register('processedFile',
                 '/home/tquast/tbMarch2018_DESY/reco/MERGED_1174.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('layerPositionFile',
                 '/afs/cern.ch/user/t/tquast/CMSSW_9_3_0/src/HGCal/CondObjects/data/layer_distances_DESY_March2018_config4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')


options.register('SkipFirstNEventsInTelescopeFile',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'SkipFirstNEventsInTelescopeFile'
                )


options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '.'
                )



options.maxEvents = -1

options.parseArguments()
print options

layerPositionFile=options.layerPositionFile

################################
process = cms.Process("TelescopeMerger")
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
                                  fileName = cms.untracked.string(options.processedFile)
)


process.daturaproducer.OutputCollectionName = cms.string("HGCalTBDATURATracks") 
process.daturaproducer.RUNDATA = cms.InputTag("source","RunData")
process.daturaproducer.inputFile = cms.string(options.TelescopeFile)
process.daturaproducer.SkipFirstNEventsInTelescopeFile = cms.int32(options.SkipFirstNEventsInTelescopeFile)
process.daturaproducer.layerPositionFile = cms.string(layerPositionFile)




process.p = cms.Path( process.daturaproducer )


process.end = cms.EndPath(process.output)
