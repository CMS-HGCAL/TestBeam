import FWCore.ParameterSet.Config as cms

process = cms.Process("unpack")
process.load('HGCal.RawToDigi.hgcaltbdigis_cfi')

process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(2), ### maybe this should be read from the file
                            fileNames=cms.untracked.vstring("file:SKIROC_RO.txt"), ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            nSpills=cms.untracked.uint32(6)
)

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")

process.plot = cms.EDAnalyzer("DigiPlotter")

process.TFileService = cms.Service("TFileService", fileName = cms.string("test_DigiPlotter_OneLayer_TB.root") )

process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string("test_Digis_OneLayer_TB.root")
        )

process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.plot)

process.outpath = cms.EndPath(process.out)
