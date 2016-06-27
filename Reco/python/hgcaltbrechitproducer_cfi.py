import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                              OutputCollectionName = cms.string(''),
                              digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_M1254.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_Check_M1254.txt'),
#				pedestalLow = cms.string(''),
#                                pedestalHigh = cms.string(''),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                              )

