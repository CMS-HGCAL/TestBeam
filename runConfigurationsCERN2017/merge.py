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

options.register('DWCFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the reconstructed DWC file.')

options.register('processedFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/rawhits/RAWHITS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('layerPositionFile',
                 'layer_distances_CERN_Hexaboard_July_6Layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('reportEvery',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '.'
                )


options.maxEvents = -1

options.parseArguments()
print options

################################
process = cms.Process("DWCMerger")
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
                                                                         'keep *_*_HGCALTBRECHITS_*',
                                                                         'keep *_*_DelayWireChambers_*',
                                                                         'keep *_*_HGCalTBDWCTracks_*',
                                                                         'keep *_*_FullRunData_*')
)


process.wirechamberproducer.OutputCollectionName = cms.string("DelayWireChambers") 
process.wirechamberproducer.RUNDATA = cms.InputTag("source","RunData")
process.wirechamberproducer.inputFile = cms.string(options.DWCFile)


process.dwctrackproducer = cms.EDProducer("DWCTrackProducer",
                                        MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                        OutputCollectionName=cms.string("HGCalTBDWCTracks"),
                                        layerPositionFile=cms.string(options.layerPositionFile)
)
process.p = cms.Path( process.wirechamberproducer*process.dwctrackproducer )


process.end = cms.EndPath(process.output)
