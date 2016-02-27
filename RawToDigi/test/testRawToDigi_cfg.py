import FWCore.ParameterSet.Config as cms

process = cms.Process("unpack")


process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(2),
                            fileNames=cms.untracked.vstring("file:SKIROC_RO.txt")

)

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")

process.hgcaldigis = cms.EDProducer("HGCalTBRawToDigi",
                                    InputLabel=cms.InputTag("source"),
                                    fedId=cms.untracked.int32(1000),
                                    electronicsMap=cms.untracked.string("HGCal/CondObjects/test/map_FNAL.txt")
                                    )

process.p =cms.Path(process.dumpRaw*process.hgcaldigis*process.dumpDigi)
