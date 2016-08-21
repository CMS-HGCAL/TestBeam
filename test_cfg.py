import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'

options.register('dataFolder',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw text input')

options.register('outputFolder',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Result of processing')

options.register('commonPrefix',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Input file to process')

options.register('runNumber',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input file to process')

options.register('chainSequence',
                 3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '1: do Digi, 3: do both Digi and Reco (only Reco not implemented so far)')

options.register('nSpills',
                 6,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of spills in run')

options.register('pedestalsHighGain',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to high gain pedestals file')

options.register('pedestalsLowGain',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to low gain pedestals file')

options.parseArguments()

if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

if not os.path.isdir(options.outputFolder):
    os.system("mkdir -p " + options.outputFolder)

process = cms.Process("unpack")
process.load('HGCal.RawToDigi.hgcaltbdigis_cfi')
process.load('HGCal.RawToDigi.hgcaltbdigisplotter_cfi')
process.load('HGCal.Reco.hgcaltbrechitproducer_cfi')
process.load('HGCal.Reco.hgcaltbrechitplotter_cfi')

process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(1), ### maybe this should be read from the file
                            #fileNames=cms.untracked.vstring("file:Raw_data_New.txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            fileNames=cms.untracked.vstring("file:%s/%s_%d.txt"%(options.dataFolder,options.commonPrefix,options.runNumber)), ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            nSpills=cms.untracked.uint32(options.nSpills)
)

process.hgcaltbdigisplotter = cms.EDAnalyzer("DigiPlotter",
                                     pedestalsHighGain=cms.untracked.string(options.pedestalsHighGain),
                                     pedestalsLowGain=cms.untracked.string(options.pedestalsLowGain)
                                     )

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")


process.output = cms.OutputModule("PoolOutputModule",
			fileName = cms.untracked.string("test_output.root")
                                 )

# process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Display.root") )
if (options.chainSequence == 1):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_%d_Digi.root"%(options.outputFolder,options.commonPrefix,options.runNumber)))
elif (options.chainSequence == 3):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_%d_Reco.root"%(options.outputFolder,options.commonPrefix,options.runNumber)))
# process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Layer.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Cluster.root") )

if (options.chainSequence == 1):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)
elif (options.chainSequence == 3):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)
# process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.FourLayerRecHitPlotterMax)
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.LayerSumAnalyzer)
process.end = cms.EndPath(process.output)
