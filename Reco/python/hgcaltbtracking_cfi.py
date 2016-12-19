import FWCore.ParameterSet.Config as cms

hgcaltbtrackanalyzer = cms.EDAnalyzer("TrackAnalyzer",
                                      HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                      HGCALTBCALOTRACKS = cms.InputTag("hgcaltbcalotracks","","unpack" ),
                                      Nlayers = cms.untracked.int32( 8 ),
                                      NSkirocsPerLayer = cms.untracked.int32( 2 ),
                                      NChannelsPerSkiroc = cms.untracked.int32( 64 ),
                                      maxChi2 = cms.untracked.double(9.48),
                                      noiseEnergyThreshold = cms.untracked.double( 9.0 ),
                                      LayerZPositions = cms.untracked.vdouble(0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4),#cern config 2
                                      #LayerZPositions = cms.untracked.vdouble(0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92),#cern config 1
                                      SensorSize = cms.untracked.int32(128)
                                      )
