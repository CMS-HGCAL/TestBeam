import FWCore.ParameterSet.Config as cms

process = cms.Process("unpack")


process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(2), ### maybe this should be read from the file
                            fileNames=cms.untracked.vstring("file:SKIROC_RO.txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED

)

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")

process.hgcaldigis = cms.EDProducer("HGCalTBRawToDigi",
                                    InputLabel=cms.InputTag("source"), ###
                                    fedId=cms.untracked.int32(1000), ## FED ID is hard coded in the HGCalTBTextSource TO BE FIXED

                                    electronicsMap=cms.untracked.string("HGCal/CondObjects/data/map_FNAL.txt"), ### electronicsMap written by hand

                                    )

process.p =cms.Path(process.dumpRaw*process.hgcaldigis*process.dumpDigi)
