import FWCore.ParameterSet.Config as cms

#hgcaltbrechitsplotter = cms.EDAnalyzer("RecHitPlotter",
#               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
#                              )
hgcaltbrechitsplotter_new = cms.EDAnalyzer("RecHitPlotter_New",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_new = cms.EDAnalyzer("RecHitPlotter_HighGain_New",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechitshighgain","","unpack" )
                              )

hgcaltbrechitsplotter_lowgain_new = cms.EDAnalyzer("RecHitPlotter_LowGain_New",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechitslowgain","","unpack" )
                              )
