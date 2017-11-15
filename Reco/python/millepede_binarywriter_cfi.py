import FWCore.ParameterSet.Config as cms


millepede_binarywriter  = cms.EDAnalyzer("MillepedeBinaryWriter",
                                binaryFile = cms.string('~/millebin.bin'),
                                nLayers = cms.int32(4),
                                useMWCReference = cms.bool(True),
                                MWCQualityCut = cms.bool(True),
                                MWCHAMBERS = cms.InputTag("source","WireChambers","unpack" ), 
                                RUNDATA = cms.InputTag("source","RunData","unpack" ), 
                                fittingMethod = cms.string("lineAnalytical")
                              )
