import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('outputFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('VariablesToPlot',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Indicate which variables are to be added to the tree. The keys here must match the ones as defined in the VariableComputation'
                )

options.register('electronicMap',
                 'map_DESY_March2018_config4_V0.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('NHexaBoards',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('hgcalLayout',
                 'layerGeom_desymarch2018_configuration4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')


options.register('layerPositionFile',
                 '/afs/cern.ch/user/t/tquast/CMSSW_9_3_0/src/HGCal/CondObjects/data/layer_distances_DESY_March2018_config4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 ''
                )

options.maxEvents = -1

options.parseArguments()
print options


electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap
hgcalLayout="HGCal/CondObjects/data/%s" % options.hgcalLayout
layerPositionFile=options.layerPositionFile


################################
process = cms.Process("analysis")
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



process.cellenergyplotting = cms.EDAnalyzer("CellEnergyPlotter",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                HGCALTBCOMMONMODENOISE = cms.InputTag("rechitproducer","HGCALTBCOMMONMODENOISEMAP" ),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                commonModeNoiseRejectionType = cms.int32(0)       #0: none, else 1-..., default: 0 
)

process.mipfindinganalysis = cms.EDAnalyzer("MIPFinder",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                MWCHAMBERS = cms.InputTag("daturaproducer","DaturaTelescopeClusters" ), 
                                DWCTRACKS = cms.InputTag("telescopetrackproducer","HGCalTBDATURATracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                HGCALTBCOMMONMODENOISE = cms.InputTag("rechitproducer","HGCALTBCOMMONMODENOISEMAP" ),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                n_bins_DWCE = cms.int32(50),
                                max_dim_x_DUT = cms.double(10),
                                max_dim_y_DUT = cms.double(10.),
                                pathsToMIPWindowFiles = cms.vstring(""),
                                commonModeNoiseRejectionType = cms.int32(0),
                                DWCs_CERNSPS = cms.untracked.bool(False) 
                              )

process.variablecomputation = cms.EDProducer("VariableComputation",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                MWCHAMBERS = cms.InputTag("daturaproducer","DaturaTelescopeClusters" ), 
                                DWCTRACKS = cms.InputTag("telescopetrackproducer","HGCalTBDATURATracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                UserRecordCollectionName=cms.untracked.string("VariableUserRecords"),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                layerPositionFile=cms.string(layerPositionFile),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                NLayers=cms.untracked.int32(options.NHexaBoards),
                                NColorsInputImage = cms.untracked.int32(-1),
                                CellEnergyCut = cms.untracked.double(0.5)
)

process.ntupelizer = cms.EDAnalyzer("NTupelizer",
                                USERRECORDS = cms.InputTag("variablecomputation","VariableUserRecords" ),
                                UserRecordKeys = cms.vstring(options.VariablesToPlot)
)

process.eventdisplay = cms.EDAnalyzer("EventDisplay",
                                RUNDATA = cms.InputTag("daturaproducer", "FullRunData" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                MWCHAMBERS = cms.InputTag("daturaproducer","DaturaTelescopeClusters" ), 
                                DWCTRACKS = cms.InputTag("telescopetrackproducer","HGCalTBDATURATracks" ), 
                                electronicsMap = cms.untracked.string(electronicMap),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                eventsToPlot=cms.vint32([1, 2, 3])
                              )



####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.cellenergyplotting*process.mipfindinganalysis*process.variablecomputation*process.ntupelizer*process.eventdisplay)

