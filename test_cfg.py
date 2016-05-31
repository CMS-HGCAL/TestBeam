import FWCore.ParameterSet.Config as cms

process = cms.Process("unpack")

#============================================================ Loading standard sequences
process.load('HGCal.StandardSequences.RawToDigi_cff')
process.load('HGCal.StandardSequences.LocalReco_cff')

process.load('HGCal.Calibration.pedestals_cfi')
process.load('HGCal.RawToDigi.hgcaltbdigisplotter_cfi')
process.load('HGCal.Reco.hgcaltbrechitplotter_cfi')

process.source = cms.Source("HGCalTBTextSource",
                            fileNames=cms.untracked.vstring("file:HGC_Output_123.txt") ### here a vector is provided
#                            fileNames=cms.untracked.vstring("file:Electron_Runs_0142076/HGC_Output_8240.txt") ### here a vector is provided
#                            fileNames=cms.untracked.vstring("file:Electron_Runs_0442076/HGC_Output_20732.txt")
#                            fileNames=cms.untracked.vstring("file:Proton_Runs_0242016/HGC_Output_8306.txt") ### here a vector is provided


)


#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_320_Reco.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_8272_Covariance.root") )

process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_123_EventDisplay.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_20798_CellHist.root") )


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
#============================================================ dumpers for debug
process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                                 dumpPayload=cms.untracked.bool(True)
                             )

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")



#============================================================ Sequences
process.debugRawSeq = cms.Sequence(process.dumpRaw)
process.DQMSeq = cms.Sequence(process.hgcaltbdigisplotter * process.hgcaltbdigisplotter_new * process.hgcaltbrechitsplotter)

#process.p =cms.Path(process.RawToDigiSeq*process.LocalRecoSeq)

#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)

#process.p =cms.Path(process.debugRawSeq * process.RawToDigiSeq ) #*  process.DQMSeq * process.pedestals )
#process.p =cms.Path(process.RawToDigiSeq * process.LocalRecoSeq * process.hgcaltbrechitsplotter_highgain_cluster ) #* process.pedestals )
#process.p =cms.Path(process.RawToDigiSeq * process.LocalRecoSeq * process.hgcaltbrechitsplotter_highgain_cluster_telescope) #* process.pedestals )
#process.p =cms.Path(process.RawToDigiSeq * process.LocalRecoSeq*process.hgcaltbrechitsplotter)
process.p =cms.Path(process.RawToDigiSeq * process.LocalRecoSeq*process.hgcaltbrechitsplotter*process.hgcaltbrechitsplotter_highgain_new) #* process.pedestals )
#process.p =cms.Path(process.RawToDigiSeq * process.LocalRecoSeq * process.hgcaltbrechitsplotter_highgain_correlation ) #* process.pedestals )
#process.p =cms.Path(process.RawToDigiSeq * process.LocalRecoSeq * process.hgcaltbrechitsplotter_highgain_correlation_cm ) #* process.pedestals )

process.end = cms.EndPath(process.output)
