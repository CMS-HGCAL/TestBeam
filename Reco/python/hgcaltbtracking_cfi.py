import FWCore.ParameterSet.Config as cms

hgcaltbtrackanalyzer = cms.EDAnalyzer("TrackAnalyzer",
                                      HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                      HGCALTBCALOTRACKS = cms.InputTag("hgcaltbcalotracks","","unpack" ),
                                      Nlayers = cms.untracked.int32( 8 ),
                                      NSkirocsPerLayer = cms.untracked.int32( 2 ),
                                      NChannelsPerSkiroc = cms.untracked.int32( 64 ),
                                      maxChi2 = cms.untracked.double(9.48),
                                      noiseEnergyThreshold = cms.untracked.double( 9.0 )
                                      )
