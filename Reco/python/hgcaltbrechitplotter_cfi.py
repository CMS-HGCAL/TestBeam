import FWCore.ParameterSet.Config as cms

hgcaltbrechitsplotter = cms.EDAnalyzer("RecHitPlotter",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )
