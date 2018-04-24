import FWCore.ParameterSet.Config as cms

daturaproducer = cms.EDProducer("HGCalTBDATURATelescopeProducer",
                            OutputCollectionName = cms.string("BeamTelescope"), 
                            RUNDATA = cms.InputTag("source","RunData","unpack"),
                            inputFile = cms.string(""),
                            SkipFirstNEventsInTelescopeFile = cms.int32(1)	#first event in the beam telescope file is a BORE event
)