import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis')
                              )
hgcaltbrechitslowgain = cms.EDProducer("HGCalTBRecHitProducerLowGain",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis')
                              )
hgcaltbrechitshighgain = cms.EDProducer("HGCalTBRecHitProducerHighGain",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis')
                              )
