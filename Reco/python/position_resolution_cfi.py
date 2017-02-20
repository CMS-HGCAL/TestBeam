import FWCore.ParameterSet.Config as cms

position_resolution_analyzer = cms.EDAnalyzer("Position_Resolution_Analyzer",
                                e_mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                                alignmentParameterFiles = cms.vstring(['']),
                                considerationMethod = cms.string('all'),
                                weightingMethod = cms.string('squaredWeighting'),
                                fittingMethod = cms.string('lineAnalytical'),
                                fitPointWeightingMethod = cms.string('none'),
                                pedestalThreshold = cms.double(2.),   
                                layers_config  = cms.int32(1),
                                ADC_per_MIP = cms.vdouble([17.31, 17.12, 16.37, 17.45, 17.31, 16.98, 16.45, 16.19, 17.55, 17.19, 16.99, 17.92, 15.95, 16.64, 16.79, 15.66]),
                                nLayers = cms.int32(8),
                                SensorSize = cms.int32(133),
                                useMWCReference = cms.bool(True),
                                RUNDATA = cms.InputTag("source","RunData","unpack" ), 
                                MWCHAMBERS = cms.InputTag("source","MultiWireChambers","unpack" ), 
                                HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                                HGCALTBCLUSTERS7 = cms.InputTag("hgcaltbclusters","7","unpack" ),
                                HGCALTBCLUSTERS19 = cms.InputTag("hgcaltbclusters","19","unpack" )
                              )
