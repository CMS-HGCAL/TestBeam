import FWCore.ParameterSet.Config as cms

hgcaltbcalotracks = cms.EDProducer("HGCalTBCaloTrackProducer",
                                   OutputCollectionName = cms.string(''),
                                   HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                                   doTrackCleaning = cms.untracked.bool( False ),
                                   maxDistanceToRecoTrack = cms.untracked.double( 1.30 ),
                                   minTouchedLayers = cms.untracked.int32( 4 )
                                 )
