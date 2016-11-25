import FWCore.ParameterSet.Config as cms

position_resolution_analyzer = cms.EDAnalyzer("Position_Resolution_Analyzer",
                                considerationMethod = cms.string('all'),
                                weightingMethod = cms.string('squaredWeighting'),
                                fittingMethod = cms.string('lineTGraphErrors'),
                                pedestalThreshold = cms.double(2.),   
                                Layer_Z_Positions = cms.vdouble([1.2, 2., 3.5, 4.3, 5.8, 6.3, 8.7, 9.5, 11.4, 12.2, 13.8, 14.6, 16.6, 17.4, 20., 20.8]),
                                ADC_per_MIP = cms.vdouble([17.24, 16.92, 17.51, 16.4, 17.35, 17.49, 16.29, 16.32]),
                                nLayers = cms.int32(8),
                                SensorSize = cms.int32(128),
                                EventsFor2DGraphs = cms.vint32([1]),
                                RUNDATA = cms.InputTag("source","RunData","unpack" ), 
                                HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )#,
                              )
