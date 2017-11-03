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

options.register('simulation',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '1=Analysis is run on simulated samples'
                 )

####################################
# Options related to the experimental setup
options.register('configuration',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '1 is the July setup (6 layers), 2 is the September setup (17 layers), 3 is the October setup (20 layers)'
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

nLayers = 6
if int(options.configuration)==1:
    nLayers = 6
elif int(options.configuration)==2:
    nLayers=17
elif int(options.configuration)==3:
    nLayers=12      


if options.simulation==1:
    rundata_tag = cms.InputTag("source", "FullRunData", "rechitproducer") 
    rechit_tag = cms.InputTag("source", "HGCALTBRECHITS", "rechitproducer")
else:
    rundata_tag = cms.InputTag("wirechamberproducer","FullRunData")
    rechit_tag = cms.InputTag("rechitproducer","HGCALTBRECHITS" )

####################################
process.energy_sum_analyzer = cms.EDAnalyzer("Energy_Sum_Analyzer",
                                ADC_per_MIP = cms.vdouble([49.3]*20*4),
                                nLayers = cms.int32(nLayers),
                                SensorSize = cms.int32(133),
                                RUNDATA = rundata_tag,
                                HGCALTBRECHITS = rechit_tag
                              )

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.energy_sum_analyzer)

