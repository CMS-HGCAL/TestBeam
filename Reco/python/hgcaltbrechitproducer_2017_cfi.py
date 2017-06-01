import FWCore.ParameterSet.Config as cms


hgcaltbrechits2017 = cms.EDProducer("HGCalTBRecHitProducer_2017",
                                OutputCollectionName = cms.string(''),
                                rawHitCollection = cms.InputTag('rawhitproducer', 'HGCALTBRAWHITS'),
                                pedestalLow = cms.string(''),
                                pedestalHigh = cms.string(''),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                                LG2HG_May2017 = cms.vdouble(10., 10., 10., 10.),   #todo: only dummy values for four skirocs are added
                                adcSaturation = cms.int32(1800),
                                mapFile_CERN_May2017 = cms.string('HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt'),
                                layers_config = cms.int32(1)        #corresponding to May 2017
                              )
