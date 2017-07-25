import FWCore.ParameterSet.Config as cms

wirechamberproducer = cms.EDProducer("HGCalTBBeamWireChamberProducer",
                            OutputCollectionName = cms.string("DelayWireChambers"), 
                            RUNDATA = cms.InputTag("source","RunData","unpack"),
                            inputFile = cms.string("")
)