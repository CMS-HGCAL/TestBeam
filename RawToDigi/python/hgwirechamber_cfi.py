import FWCore.ParameterSet.Config as cms

wirechamber_producer = cms.EDProducer("HGCalTBBeamWireChamberProducer",
                            MWCHAMBERS = cms.InputTag("source","WireChambers","unpack"),
                            RUNDATA = cms.InputTag("source","RunData","unpack"),
                            inputFile = cms.string("")
)