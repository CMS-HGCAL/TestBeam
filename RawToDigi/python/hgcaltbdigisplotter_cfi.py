import FWCore.ParameterSet.Config as cms

hgcaltbdigisplotter = cms.EDAnalyzer("DigiPlotter",
                                     pedestalsHighGain=cms.untracked.string(""),
                                     pedestalsLowGain=cms.untracked.string("")
                                     )
hgcaltbdigisplotter_new = cms.EDAnalyzer("DigiPlotter_New")
hgcaltbdigisplotter_ped = cms.EDAnalyzer("DigiPlotter_Ped")
