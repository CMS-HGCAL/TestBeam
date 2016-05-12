import FWCore.ParameterSet.Config as cms

hgcaltbtracks = cms.EDProducer("FNALTelescopeRawToTracks",
                            InputLabel = cms.InputTag("source"), ###
                            fedIds = cms.untracked.vint32(99), ## list of FEDs should be know from the setup
)

