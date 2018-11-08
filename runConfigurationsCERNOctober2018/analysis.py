import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'file containing input')

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
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('NHexaBoards',
                -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('NLayers',
                -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of layers for analysis.'
                )

options.register('hgcalLayout',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')


options.register('layerPositionFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('isSimulation',
                0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 ''
                )

options.register('runDWCCorrelation',
                0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 ''
                )

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

if options.isSimulation==0:
    rundata_tag = cms.InputTag("wirechamberproducer", "FullRunData" )
    rechit_tag = cms.InputTag("rechitproducer","HGCALTBRECHITS" )
    commonModeNoise_tag = cms.InputTag("rechitproducer","HGCALTBCOMMONMODENOISEMAP" )
    dwc_tag = cms.InputTag("wirechamberproducer","DelayWireChambers" )
    dwc_track_tag = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" )
else:
    rundata_tag = cms.InputTag("source", "FullRunData" )
    rechit_tag = cms.InputTag("source","HGCALTBRECHITS" )
    commonModeNoise_tag = cms.InputTag("","" )
    dwc_tag = cms.InputTag("source","DelayWireChambers" )
    dwc_track_tag = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" )



process.variablecomputation = cms.EDProducer("VariableComputation",
                                RUNDATA = rundata_tag,  
                                MWCHAMBERS = dwc_tag, 
                                DWCTRACKS = dwc_track_tag,                                 
                                HGCALTBRECHITS = rechit_tag,
                                UserRecordCollectionName=cms.untracked.string("VariableUserRecords"),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                layerPositionFile=cms.string(layerPositionFile),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                NLayers=cms.untracked.int32(options.NLayers),
                                NColorsInputImage = cms.untracked.int32(-1),
                                CellEnergyCut = cms.untracked.double(0.5)
)

process.observablentupler = cms.EDAnalyzer("NTupelizer",
                                USERRECORDS = cms.InputTag("variablecomputation","VariableUserRecords" ),
                                UserRecordKeys = cms.vstring(options.VariablesToPlot)
)

process_chain = process.variablecomputation*process.observablentupler

if options.runDWCCorrelation!=0:
    process.correlationanalysis = cms.EDAnalyzer("DWCCorrelator",
                                    RUNDATA = rundata_tag, 
                                    MWCHAMBERS = dwc_tag, 
                                    DWCTRACKS = dwc_track_tag, 
                                    HGCALTBRECHITS = rechit_tag,
                                    HGCALTBCOMMONMODENOISE = commonModeNoise_tag,
                                    ElectronicMap = cms.untracked.string(electronicMap),
                                    DetectorLayout=cms.untracked.string(hgcalLayout),
                                    NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                    n_bins_DWCE = cms.int32(50),
                                    max_dim_x_DUT = cms.double(5),
                                    max_dim_y_DUT = cms.double(5),
                                    pathsToMIPWindowFiles = cms.vstring([""]),
                                    commonModeNoiseRejectionType = cms.int32(0)       #0: none, else 1-..., default: 0 
                                  )
    process_chain = process_chain * process.correlationanalysis


process.p = cms.Path(process_chain)

