import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/home/tquast/tb2017/reconstructedFiles/reco_1511.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('outputFile',
                 '/home/tquast/tb2017/analysis/MIPs_1259.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where analysis output are stored')


options.register('electronicMap',
                 'map_CERN_Hexaboard_September_17Sensors_7EELayers_10FHLayers_V1.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('NHexaBoards',
                10,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('NLayers',
                10,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of layers of HGCal prototype for analysis.'
                )

options.register('NDWCs',
                4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of Delay Wire chambers in beam line.'
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

################################
process = cms.Process("mipfindinganalysis")
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

process.numpyconverter = cms.EDAnalyzer("NumpyConverter",
                                RUNDATA = cms.InputTag("wirechamberproducer", "FullRunData" ), 
                                MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                DWCTRACKS = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                electronicMap = cms.untracked.string(electronicMap),
                                NHexaBoards=cms.untracked.uint32(options.NHexaBoards),
                                NLayers=cms.untracked.uint32(options.NLayers),
                                NDWCs=cms.untracked.uint32(options.NDWCs),
                                Sensorsize=cms.untracked.uint32(128),
                                outputFilePath=cms.string(options.outputFile),
                                eventIdentifier=cms.untracked.string("event"),
                                rechitsIdentifier=cms.untracked.string("rechits"),
                                dwcsIdentifier=cms.untracked.string("dwcs"),
                                dwcReferenceIdentifier=cms.untracked.string("dwcReference")
                              )


####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.numpyconverter)

