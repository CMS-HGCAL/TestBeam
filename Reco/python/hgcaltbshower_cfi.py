import FWCore.ParameterSet.Config as cms

hgcaltbshower = cms.EDAnalyzer("ShowerAnalyzer",
                               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                               HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                               HGCALTBTRACKS = cms.InputTag("hgcaltbcalotracks","","unpack" ),
                               Nlayers = cms.untracked.int32( 8 ),
                               NSkirocsPerLayer = cms.untracked.int32( 2 ),
                               SensorSize = cms.untracked.int32( 128 ),
                               CERN_8layers_config = cms.untracked.int32( 0 ),
                               maxTransverseProfile = cms.untracked.double( 250 )
                               )
