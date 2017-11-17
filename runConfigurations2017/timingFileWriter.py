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

options.register('timingFile',
                 '/eos/user/t/tquast/data/Testbeam/July2017/Timing/HexaData_Run1303_TIMING_RDOUT_ORM0.txt.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


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
process = cms.Process("timingfilewriter")
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


process.timingfilewriter = cms.EDAnalyzer("HGCalTBTimingFileWriter",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        TimingFilePath=cms.untracked.string(options.timingFile)
)



process.p = cms.Path( process.timingfilewriter )


