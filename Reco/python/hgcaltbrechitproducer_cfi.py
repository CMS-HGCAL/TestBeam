import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_Test_2Layer.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_Test_2Layer.txt'),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                              )

