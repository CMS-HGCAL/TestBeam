import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/home/tquast/tb2017/reconstructedFiles/reco_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('outputFile',
                 '/home/tquast/tb2017/analysis/energyReco_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )


####################################
# Options related to the experimental setup
options.register('configuration',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 '1 is the July setup'
                 )



options.maxEvents = -1


options.parseArguments()
print options



################################
process = cms.Process("positionresolutionanalysis")
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

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))

####################################
process.energy_sum_analyzer = cms.EDAnalyzer("Energy_Sum_Analyzer",
                                layers_config  = cms.int32(1),
                                ADC_per_MIP = cms.vdouble([1.]*4*10),
                                nLayers = cms.int32(6),
                                SensorSize = cms.int32(133),
                                RUNDATA = cms.InputTag("wirechamberproducer","FullRunData" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                              )


####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.energy_sum_analyzer)

