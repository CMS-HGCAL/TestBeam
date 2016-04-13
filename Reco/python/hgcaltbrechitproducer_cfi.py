import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string(''),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_8272_Mean.txt'),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                              )

