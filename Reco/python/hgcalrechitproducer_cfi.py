import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis')
                              )

