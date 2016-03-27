import FWCore.ParameterSet.Config as cms

hgcaltbdigis = cms.EDProducer("HGCalTBRawToDigi",
                            InputLabel=cms.InputTag("source"), ###
                            fedIds=cms.untracked.vint32(85, 88), ## FED ID is hard coded in the HGCalTBTextSource TO BE FIXED
                            
                            electronicsMap=cms.untracked.string("HGCal/CondObjects/data/map_FNAL_1.txt"), ### electronicsMap written by hand
                            
)

