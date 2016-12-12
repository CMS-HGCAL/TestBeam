import FWCore.ParameterSet.Config as cms

hgcaltbcalotracks = cms.EDProducer("HGCalTBCaloTrackProducer",
                                   OutputCollectionName = cms.string(''),
                                   HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                                   HGCALTBRECHITS = cms.InputTag('hgcaltbrechits',"","unpack"),
                                   LayerZPositions= cms.untracked.vdouble(0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4),#cern config 2
                                   doTrackCleaning = cms.untracked.bool( False ),
                                   maxDistanceToRecoTrack = cms.untracked.double( 2.0 ),
                                   minTouchedLayers = cms.untracked.int32( 4 ),
                                   minEnergy = cms.untracked.double( 9.0 ),
                                   maxEnergy = cms.untracked.double( 60.0 )
                                 )
