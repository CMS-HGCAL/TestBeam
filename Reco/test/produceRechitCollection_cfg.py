import FWCore.ParameterSet.Config as cms
processName = "HGC"
process = cms.Process(processName)

process.load('FWCore.MessageService.MessageLogger_cfi')


process.source = cms.Source ("EmptySource")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard')
options.maxEvents = 200
process.MessageLogger.cerr.FwkReport.reportEvery = 1 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)



process.RechitCollection = cms.EDProducer("produceRechitCollection",
                OutputCollectionName1 = cms.string("HGCRechitCollection")
)

process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string("test_RecHits_OneLayer_TB.root")
        )


# Defines which modules and sequences to run
process.mypath = cms.Path(process.RechitCollection)

#process.schedule = cms.Schedule(process.mypath)

# A list of analyzers or output modules to be run after all paths have been run.
process.outpath = cms.EndPath(process.out)
