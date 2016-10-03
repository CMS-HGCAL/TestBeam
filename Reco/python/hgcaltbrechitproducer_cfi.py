import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                                OutputCollectionName = cms.string(''),
                                digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_L8.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_L8.txt'),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                                LG2HG = cms.double(9.),
                                gainThr = cms.double(1800.)
                                )
