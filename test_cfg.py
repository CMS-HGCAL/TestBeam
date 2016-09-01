import FWCore.ParameterSet.Config as cms

process = cms.Process("unpack")
process.load('HGCal.RawToDigi.hgcaltbdigis_cfi')
process.load('HGCal.RawToDigi.hgcaltbdigisplotter_cfi')
process.load('HGCal.Reco.hgcaltbrechitproducer_cfi')
process.load('HGCal.Reco.hgcaltbrechitplotter_cfi')

RUNNUMBER  = "843"
process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(1),####provide file name below
                            fileNames=cms.untracked.vstring("file:/afs/cern.ch/work/r/rslu/public/HGC_TB_data_Sep2016/HGCRun_Output_"+RUNNUMBER+".txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
#                            fileNames=cms.untracked.vstring("file:/afs/cern.ch/work/r/rslu/public/HGC_TB_data_Sep2016/PED_Output_000"+RUNNUMBER+".txt")
)

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")


process.output = cms.OutputModule("PoolOutputModule",
			fileName = cms.untracked.string("test_output.root")
                                 )

#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_Cluster_"+RUNNUMBER+".root")) ### Analyzed output file with histograms
process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_Digi_"+RUNNUMBER+".root")) ### Analyzed output file with histograms
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_Reco_"+RUNNUMBER+".root")) ### Analyzed output file with histograms
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_EventDisplay_"+RUNNUMBER+".root")) ### Analyzed output file with histograms

########Activate this to produce event displays#########################################
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)

################Not needed for DQM purposes, produces digi histograms for each channel, and the pedestal txt file needed for Digi->Reco
process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)

################Produces Reco histograms for each channel as well as a scatter plot of the Reco per channel#############
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)

#################Produces Clusters of Recos(7cells, 19cells and all cells(full hexagons only))################
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.LayerSumAnalyzer)


process.end = cms.EndPath(process.output)
