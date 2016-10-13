import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                                OutputCollectionName = cms.string(''),
                                digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_L8.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_L8.txt'),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                                LG2HG = cms.vdouble(10.2, 10.2, 10., 10., 9.8, 8.8, 9.7, 9.7, 9.7, 9.7, 9.8, 9.8, 9.8, 9.8, 9.2, 9.2),
                                adcSaturation = cms.int32(1800)
                              )
