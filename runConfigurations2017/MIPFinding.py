import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/home/tquast/tb2017/reconstructedFiles/reco_1259.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('outputFile',
                 '/home/tquast/tb2017/analysis/MIPs_1259.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where analysis output are stored')

options.register('pathsToMIPWindowFiles',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the file containing the corresponding cuts on the reconstructed impact positions on DWC E.'
                )

options.register('electronicMap',
                 'map_CERN_Hexaboard_July_6Layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('hgcalLayout',
                 'layerGeom_oct2017_h2_17layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')

options.register('NHexaBoards',
                10,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('NBins',
                60,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Bins for the DWC search window.'
                )

options.register('DimXDUT',
                30.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'MIP search window on DWCE in x-coordinate [mm].'
                )

options.register('DimYDUT',
                30.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'MIP search window on DWCE in y-coordinate [mm].'
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
hgcalLayout="HGCal/CondObjects/data/%s" % options.hgcalLayout

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


process.mipfindinganalysis = cms.EDAnalyzer("DWCCorrelator",
                                RUNDATA = cms.InputTag("wirechamberproducer", "FullRunData" ), 
                                MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                DWCTRACKS = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                HGCALTBCOMMONMODENOISE = cms.InputTag("rechitproducer","HGCALTBCOMMONMODENOISEMAP" ),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                n_bins_DWCE = cms.int32(options.NBins),
                                max_dim_x_DUT = cms.double(options.DimXDUT),
                                max_dim_y_DUT = cms.double(options.DimYDUT),
                                pathsToMIPWindowFiles = cms.vstring(options.pathsToMIPWindowFiles),
                                commonModeNoiseRejectionType = cms.int32(0)       #0: none, else 1-..., default: 0 

                              )


####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.mipfindinganalysis)

