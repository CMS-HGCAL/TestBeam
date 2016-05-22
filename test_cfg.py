import FWCore.ParameterSet.Config as cms

process = cms.Process("unpack")

#============================================================ Loading standard sequences
process.load('HGCal.StandardSequences.RawToDigi_cff')
process.load('HGCal.StandardSequences.LocalReco_cff')

process.load('HGCal.Calibration.pedestals_cfi')
process.load('HGCal.RawToDigi.hgcaltbdigisplotter_cfi')
process.load('HGCal.Reco.hgcaltbrechitplotter_cfi')

process.source = cms.Source("HGCalTBTextSource",
                            fileNames=cms.untracked.vstring("NOFILE"), ### here a vector is provided, NOFILE should be given otherwise it would not run at all
                            hgcalData = cms.untracked.PSet(
        fileNames=cms.untracked.vstring("file:HGC_Output_112.txt"), ### here a vector is provided
                  ),
                            telescopeData = cms.untracked.PSet(
        fileNames = cms.untracked.vstring("NOFILE"), 
        ),
                            )


process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Telescope_20532_751.root") )

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
process.DQMSeq = cms.Sequence(process.hgcaltbdigisplotter * process.hgcaltbrechitsplotter)

#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.hgcaltbdigisplotter*process.LocalRecoSeq r)

#process.p =cms.Path(process.dumpRaw*process.hgcaltbdigis*process.dumpDigi*process.hgcaltbdigisplotter)

#process.p =cms.Path(process.debugRawSeq * process.RawToDigiSeq ) #*  process.DQMSeq * process.pedestals )
#process.p =cms.Path(process.RawToDigiSeq*process.LocalRecoSeq*process.hgcaltbrechitsplotter_highgain_new) # * process.LocalRecoSeq *  process.DQMSeq ) #* process.pedestals )
process.p =cms.Path(process.RawToDigiSeq*process.LocalRecoSeq*process.hgcaltbdigisplotter_new) # * process.LocalRecoSeq *  process.DQMSeq ) #* process.pedestals )
process.end = cms.EndPath(process.output)

