import FWCore.ParameterSet.Config as cms

imaging_tuple_writer = cms.EDAnalyzer("Imaging_Tuple_Writer",
                                RUNDATA = cms.InputTag("source","RunData","unpack" ), 
                                HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                              )
