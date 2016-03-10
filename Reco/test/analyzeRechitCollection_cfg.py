import FWCore.ParameterSet.Config as cms

process = cms.Process("Plot")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/test_RecHits_OneLayer_TB.root'
                )
                            )

process.plot = cms.EDAnalyzer("RecHitPlotter",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

process.TFileService = cms.Service("TFileService", fileName = cms.string("test_RecHitPlotter_OneLayer_TB.root") )

process.p = cms.Path(process.plot)
