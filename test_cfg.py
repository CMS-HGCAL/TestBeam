import FWCore.ParameterSet.Config as cms

process = cms.Process("unpack")
process.load('HGCal.RawToDigi.hgcaltbdigis_cfi')
process.load('HGCal.RawToDigi.hgcaltbdigisplotter_cfi')
process.load('HGCal.Reco.hgcaltbrechitproducer_cfi')
process.load('HGCal.Reco.hgcaltbrechitplotter_cfi')

process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(1), ### maybe this should be read from the file
#                          fileNames=cms.untracked.vstring("file:Proton_Runs_0342016/HGC_Output_8336.txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                          fileNames=cms.untracked.vstring("file:Proton_Runs_0242016/HGC_Output_8272.txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
)

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")

process.output = cms.OutputModule("PoolOutputModule",
                                  compressionAlgorithm = cms.untracked.string('LZMA'),
                                  compressionLevel = cms.untracked.int32(4),
                                  dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('USER'),
        filterName = cms.untracked.string('')
        ),
                                  #dropMetaData = cms.untracked.string('ALL'),
                                  eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
                                 fastCloning = cms.untracked.bool(False),
                                 fileName = cms.untracked.string('test_output.root'), #options.output),
#                                 outputCommands = process.MICROAODSIMEventContent.outputCommands,
#                                 overrideInputFileSplitLevels = cms.untracked.bool(True),
#                                 SelectEvents = SelectEventsPSet
                                 )

#process.TFileService = cms.Service("TFileService", fileName = cms.string("test_DigiAndRechitPlotter_TB_8308_RecHits_ADClt20_FullHex_Xgt0_Ygt0.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("test_Pedestal_8272_Correlation.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("test_8336_HighGain.root") )


#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.hgcaltbrechits*process.hgcaltbrechitsplotter*process.hgcaltbrechitsplotter_new)
process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechitshighgain*process.hgcaltbrechitsplotter_highgain_correlation)
#process.p =cms.Path(process.hgcaltbdigis*process.dumpDigi)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.hgcaltbrechitslowgain*process.hgcaltbrechitsplotter_lowgain_new)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.hgcaltbrechits*process.hgcaltbrechitsplotter*process.hgcaltbrechitsplotter_lowgain_new)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.hgcaltbrechits*process.hgcaltbrechitsplotter)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.hgcaltbdigisplotter)
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.hgcaltbdigisplotter_new)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.hgcaltbdigisplotter_highgain_new)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.hgcaltbdigisplotter_ped_profile)
#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.hgcaltbdigisplotter_HighGain_ped_profile);


process.end = cms.EndPath(process.output)
