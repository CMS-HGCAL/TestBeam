import FWCore.ParameterSet.Config as cms

hgcaltbdigis = cms.EDProducer("HGCalTBRawToDigi",
                            InputLabel=cms.InputTag("source"), ###
                            fedIds=cms.untracked.vint32(0, 1), ## FED ID is hard coded in the HGCalTBTextSource TO BE FIXED
                            
                            electronicsMap=cms.untracked.string("HGCal/CondObjects/data/map_FNAL_2.txt"), ### electronicsMap written by hand
                            
)

