import FWCore.ParameterSet.Config as cms


dwc_ntupelizer  = cms.EDAnalyzer("DWC_NTupelizer",
								MWCHAMBERS = cms.InputTag("source","WireChambers","unpack" ),
                                RUNDATA = cms.InputTag("source","RunData","unpack" ),
                                writeMinimal = cms.bool(False)
                              )
