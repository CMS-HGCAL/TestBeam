import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_Test_1Layer.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_Test_1Layer.txt'),  
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                              )

hgcaltbsimrechits = cms.EDProducer("HGCalTBSimRecHitProducer",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag("source","","HGC"),
                                pedestalLow = cms.string(''),
                                pedestalHigh = cms.string(''),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                              )

