import FWCore.ParameterSet.Config as cms

millepede_binarywriter  = cms.EDAnalyzer("MillepedeBinaryWriter",
                                binaryFile = cms.string('~/millebin.bin'),
                                considerationMethod = cms.string('all'),
                                weightingMethod = cms.string('squaredWeighting'),
                                fittingMethod = cms.string('lineAnalytical'),
                                fitPointWeightingMethod = cms.string('none'),
                                pedestalThreshold = cms.double(2.),   
                                layers_config  = cms.int32(1),
                                ADC_per_MIP = cms.vdouble([17.24, 16.92, 17.51, 16.4, 17.35, 17.49, 16.29, 16.32]),
                                nLayers = cms.int32(8),
                                SensorSize = cms.int32(128),
                                totalEnergyThreshold = cms.double(1000.),
                                RUNDATA = cms.InputTag("source","RunData","unpack" ), 
                                HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                                HGCALTBCLUSTERS7 = cms.InputTag("hgcaltbclusters","7","unpack" ),
                                HGCALTBCLUSTERS19 = cms.InputTag("hgcaltbclusters","19","unpack" )
                              )
