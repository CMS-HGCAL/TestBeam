import FWCore.ParameterSet.Config as cms

hgcaltbdigis = cms.EDProducer("HGCalTBRawToDigi",
                            InputLabel=cms.InputTag("source"), ###
                            #InputLabel=cms.InputTag("test"), ###
                            fedId=cms.untracked.int32(1000), ## list of FEDs should be know from the setup
                            
                            electronicsMap=cms.untracked.string("HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt"), ### electronicsMap written by hand
                            
)

BadSpillFilter = cms.EDFilter("Bad_Spill_Filter",
                              layers_config = cms.int32(-1),
                              nameCFG1 = cms.string('CondObjects/data/Bad_Run_Spill_CFG1.txt'),
                              nameCFG2 = cms.string('CondObjects/data/Bad_Run_Spill_CFG2.txt')

                              )
