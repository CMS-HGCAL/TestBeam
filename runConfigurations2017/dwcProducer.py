import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/rechits/RECHITS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('DWCFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/reconstructed_DWC_data_full.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the reconstructed DWC file.')

options.register('processedFile',
                 '/home/tquast/tb2017/reconstructedFiles/reco_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('reportEvery',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )


options.maxEvents = -1


options.parseArguments()
print options



################################
process = cms.Process("dwcproducer")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s"%options.dataFile)
)

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################



process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.processedFile),                                  
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*',
                                                                         'keep *_*_DelayWireChambers_*',
                                                                         'keep *_*_FullRunData_*')
)

#Wire chamber producer
process.wirechamberproducer.OutputCollectionName = cms.string("DelayWireChambers") 
process.wirechamberproducer.RUNDATA = cms.InputTag("source","RunData")
process.wirechamberproducer.inputFile = cms.string(options.DWCFile)


####################################
#add skip event exception which might occur for simulated samples because the last event is not properly passed forward
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)



process.p = cms.Path( process.wirechamberproducer)
process.end = cms.EndPath(process.output)

