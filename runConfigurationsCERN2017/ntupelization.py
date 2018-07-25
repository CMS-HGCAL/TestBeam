import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/home/tquast/tbJune2018_H2/reco/merged_000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'file containing input')

options.register('outputFile',
                 '/home/tquast/tbJune2018_H2/ntuples/ntuple_000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('electronicMap',
                 "map_CERN_Hexaboard_June_28Sensors_28EELayers_V0.txt",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('NHexaBoards',
                28,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('hgcalLayout',
                 'layerGeom_june2018_h2_28layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')

options.register('layerPositionFile',
                 '/afs/cern.ch/user/t/tquast/CMSSW_9_3_0/src/HGCal/CondObjects/data/layer_distances_CERN_Hexaboard_June2018_28Layers_dummy.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('isSimulation',
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
process = cms.Process("ntuples")
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
    dwc_tag = cms.InputTag("wirechamberproducer","DelayWireChambers" )
    dwc_track_tag = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" )
else:
    rundata_tag = cms.InputTag("source", "FullRunData" )
    rechit_tag = cms.InputTag("source","HGCALTBRECHITS" )
    dwc_tag = cms.InputTag("source","DelayWireChambers" )
    dwc_track_tag = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" )



process.eventdisplay = cms.EDAnalyzer("EventDisplay",
                                RUNDATA = rundata_tag, 
                                HGCALTBRECHITS = rechit_tag,
                                electronicsMap = cms.untracked.string(electronicMap),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                eventsToPlot=cms.vint32(range(1, 4))
                              )


process.rechitntupler = cms.EDAnalyzer("RecHitNtupler",
                                       InputCollection=rechit_tag,
                                       RUNDATA = rundata_tag,
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       layerPositionFile = cms.untracked.string(layerPositionFile),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(True),
                                       MipThreshold=cms.untracked.double(5.0),
                                       NoiseThreshold=cms.untracked.double(0.0)
)

process.trackimpactntupler = cms.EDAnalyzer("ImpactPointNtupler",
                                       extrapolationDevice=cms.untracked.string("DWC"),
                                       DWCTrackToken = dwc_track_tag,
                                       DATURATelescopeData = cms.InputTag("","" ),
                                       RUNDATA = rundata_tag,
                                       nLayers=cms.untracked.int32(options.NHexaBoards),
)


####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.eventdisplay*process.rechitntupler*process.trackimpactntupler)

