import FWCore.ParameterSet.Config as cms

hgcaltbdigis = cms.EDProducer("HGCalTBRawToDigi",
                            InputLabel=cms.InputTag("source"), ###
                            fedId=cms.untracked.int32(1000), ## list of FEDs should be know from the setup
                            
                            electronicsMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt"), ### electronicsMap written by hand
                            
)

