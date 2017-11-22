import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 106,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('outputFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('electronicMap',
                 'HGCal/CondObjects/data/map_CERN_Hexaboard_September_17Sensors_7EELayers_10FHLayers_V1.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')
# Available options: 
# "HGCal/CondObjects/data/map_CERN_Hexaboard_July_6Layers.txt"
# "HGCal/CondObjects/data/map_CERN_Hexaboard_September_17Sensors_7EELayers_10FHLayers_V1.txt" # end of september
# "HGCal/CondObjects/data/map_CERN_Hexaboard_October_17Sensors_5EELayers_6FHLayers_V1.txt" # october 18-22, 1st conf
# "HGCal/CondObjects/data/map_CERN_Hexaboard_October_20Sensors_5EELayers_7FHLayers_V1.txt" # october 18-22, 2nd conf

options.register('hgcalLayout',
                 'HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to hgcal layout file')
# Available options: 
# "HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt"
# "HGCal/CondObjects/data/layerGeom_oct2017_h6_17layers.txt"
# "HGCal/CondObjects/data/layerGeom_oct2017_h6_20layers.txt"

options.maxEvents = -1
options.output = "cmsswEvents_Cluster.root"

options.parseArguments()
print options
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

################################
process = cms.Process("rechitanalysis")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s/cmsswEvents_Run%d_RecHit.root"%(options.dataFolder,options.runNumber))
)

filename = options.outputFolder+"/HexaOutput_"+str(options.runNumber)+".root"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(filename)
)
process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*',
                                                                         'keep *_*_HGCALTBCLUSTERS_*',
                                                                         'keep *_*_HGCALTBCLUSTERS7_*',
                                                                         'keep *_*_HGCALTBCLUSTERS19_*',
                                  )
)

process.clusterproducer = cms.EDProducer("HGCalTBClusterProducer",
                                        OutputCollectionName = cms.string('HGCALTBCLUSTERS'),
                                        InputCollection = cms.InputTag('rechitproducer','HGCALTBRECHITS'),
                                        ElectronicsMap = cms.untracked.string(options.electronicMap),
                                        DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                        SensorSize=cms.untracked.int32(128),
                                        runDynamicCluster = cms.untracked.bool(True),
                                        runCluster7 = cms.untracked.bool(True),
                                        runCluster19 = cms.untracked.bool(True),
                                        minEnergy = cms.untracked.double(2.0),

)

process.rechitplotter = cms.EDAnalyzer("RecHitPlotter",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       ElectronicMap=cms.untracked.string(options.electronicMap),
                                       DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(True),
                                       MipThreshold=cms.untracked.double(5.0),
                                       NoiseThreshold=cms.untracked.double(2.0)
)

process.showeranalyzer = cms.EDAnalyzer("ShowerAnalyzer",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       InputClusterCollection=cms.InputTag("clusterproducer","HGCALTBCLUSTERS"),
                                       InputCluster7Collection=cms.InputTag("clusterproducer","HGCALTBCLUSTERS7"),
                                       InputCluster19Collection=cms.InputTag("clusterproducer","HGCALTBCLUSTERS19"),
                                       DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       NoiseThreshold=cms.untracked.double(2.0)
)


process.p = cms.Path( process.clusterproducer*process.showeranalyzer )

process.end = cms.EndPath(process.output)
