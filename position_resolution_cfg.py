import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys


options = VarParsing.VarParsing('standard')

####################################
# Options for reading in the data
options.register('fileDirectory',
                '/eos/user/t/tquast/data/Testbeam/July2017/DWC/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to the file from which the DWCs are read.'
                )

options.register('reportEvery',
                50000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.parseArguments()

#files = ["dwc_run_10.root", "dwc_run_12.root", "dwc_run_13.root", "dwc_run_14.root", "dwc_run_15.root"]
runs = [993, 994, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1012, 1013, 1025, 1027, 1028, 1029, 1030, 1032, 1033, 1035, 1036, 1039, 1045]
files = ["dwc_run_%s.root" % run for run in runs]


################################
# Setting an upper limit for the events to be processed, e.g. for debugging
options.maxEvents = -1
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
####################################


####################################
# Initialize the data read-in plugins
process.source = cms.Source("HGCalTBWireChamberSource",
    OutputCollectionName = cms.string("WireChambers"), 
    fileNames = cms.untracked.vstring(["file:%s/%s"%(options.fileDirectory, file) for file in files]),
    performAlignment = cms.untracked.bool(False),
    alignmentParamaterFile = cms.untracked.string("/tmp/millepede.res") 

)


process.millepede_binarywriter.binaryFile = cms.string('/afs/cern.ch/user/t/tquast/millepede.bin')
process.millepede_binarywriter.nLayers = cms.int32(4)
process.millepede_binarywriter.MWCQualityCut = cms.bool(True)
process.millepede_binarywriter.makeTree = cms.untracked.bool(True)
process.millepede_binarywriter.MWCHAMBERS = cms.InputTag("source","WireChambers","unpack" )
process.millepede_binarywriter.fittingMethod = cms.string("lineAnalytical")
process.millepede_binarywriter.binaryFile = cms.string("/tmp/millepede.bin")
                              


#tree file:
process.TFileService = cms.Service("TFileService", fileName = cms.string("outfile_DWCs.root"))

####################################
#add skip event exception which might occur for simulated samples because the last event is not properly passed forward
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.p = cms.Path(process.millepede_binarywriter)

