import FWCore.ParameterSet.Config as cms

hgcaltbclusters = cms.EDProducer("HGCalTBClusterProducer",
                                 OutputCollectionName = cms.string(''),
                                 rechitCollection = cms.InputTag('hgcaltbrechits',"","unpack"),
                                 LayerZPositions= cms.untracked.vdouble(0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4),#cern config 2
                                 #LayerZPositions= cms.untracked.vdouble(0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92),#cern config 1
                                 reversedLayers = cms.untracked.vint32(5,7),
                                 minEnergy = cms.untracked.double(0)
                                 )
