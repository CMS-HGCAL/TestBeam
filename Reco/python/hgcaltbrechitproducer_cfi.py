import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                              OutputCollectionName = cms.string(''),
                             digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_L16_Check.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_L16_Check.txt'),
#				pedestalLow = cms.string(''),
#                                pedestalHigh = cms.string(''),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                                adcSaturation = cms.uint32(2000),
                              )

