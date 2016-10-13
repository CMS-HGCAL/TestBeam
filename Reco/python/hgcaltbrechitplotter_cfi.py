import FWCore.ParameterSet.Config as cms

hgcaltbrechitsplotter = cms.EDAnalyzer("RecHitPlotter",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_new = cms.EDAnalyzer("RecHitPlotter_HighGain_New",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
               HGCALTBTRACKS  = cms.InputTag("hgcaltbtracks","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_cluster = cms.EDAnalyzer("RecHitPlotter_HighGain_Cluster",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_cluster_telescope = cms.EDAnalyzer("RecHitPlotter_HighGain_Cluster_Telescope",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

hgcaltbsimrechitsplotter_highgain_new = cms.EDAnalyzer("SimRecHitPlotter_HighGain_New",
               HGCALTBRECHITS = cms.InputTag("hgcaltbsimrechits","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_correlation = cms.EDAnalyzer("RecHitPlotter_HighGain_Correlation",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )
hgcaltbrechitsplotter_highgain_correlation_cm = cms.EDAnalyzer("RecHitPlotter_HighGain_Correlation_CM",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

FourLayerRecHitPlotterMax = cms.EDAnalyzer("FourLayerRecHitPlotterMax",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

LayerSumAnalyzer = cms.EDAnalyzer("Layer_Sum_Analyzer",
                                  HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                  CERN_8layers_config = cms.int32(0)
                              )
