import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'

options.register('dataFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw text input')

options.register('outputFolder',
                 '/tmp/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Result of processing')

options.register('runNumber',
                 850,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input file to process')

options.register('runType',
                 'Unknown',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Type of run: HGCRun for run with beam on, PED for pedestal run, Unknown otherwise')

options.register('chainSequence',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '0: if runType is PED then do Digi, if runType is HGC_Run then do Digi and Reco (not implemented yet); 1: do Digi; 2: only Reco (not implemented yet); 3: Digi + highgain_correlation_cm; 4: event display sequence; 5: highgain_correlation_cm + event display sequence')

options.register('nSpills',
                 15,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of spills in run')

options.register('pedestalsHighGain',
                 'CondObjects/data/pedLowGain1024.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to high gain pedestals file')

options.register('pedestalsLowGain',
                 'CondObjects/data/pedHighGain1024.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to low gain pedestals file')

options.register('configuration',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 '1 if 8Layers with 5X0 sampling the center of the shower only; 2 if 8Layers with 25X0 sampling up to the tail of the shower')



options.output = "test_output.root"
options.maxEvents = -1

options.parseArguments()

if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

if not os.path.isdir(options.outputFolder):
    os.system("mkdir -p " + options.outputFolder)

if (options.runType != "PED" and options.runType != "HGCRun"):
    sys.exit("Error: only runtypes PED and HGCRun supported for now; given runType was %s"%(options.runType))

if (options.runType == "PED"):
    if (os.path.isfile(options.pedestalsHighGain) or os.path.isfile(options.pedestalsLowGain)):
        sys.exit("Error: Run %d is a pedestals run. The arguments pedestalsHighGain = %s and pedestalsLowGain = %s should be paths that do not lead to an existing file."%(options.runNumber, options.pedestalsHighGain, options.pedestalsLowGain))
            


################################
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################
process.load('HGCal.StandardSequences.RawToDigi_cff')
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.dqm_cff')
process.load('HGCal.StandardSequences.TrackingReco_cff')
process.load('HGCal.StandardSequences.ShowerReco_cff')

process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(options.runNumber),
                            fileNames=cms.untracked.vstring("file:%s/%s_Output_%06d.txt"%(options.dataFolder,options.runType,options.runNumber)),
                            nSpills=cms.untracked.uint32(options.nSpills),
)

######
process.hgcaltbdigis.electronicsMap = cms.untracked.string("HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt")

process.hgcaltbdigisplotter.pedestalsHighGain = cms.untracked.string(options.pedestalsHighGain)
process.hgcaltbdigisplotter.pedestalsLowGain  = cms.untracked.string(options.pedestalsLowGain)

process.hgcaltbrechits.pedestalLow = cms.string(options.pedestalsLowGain)
process.hgcaltbrechits.pedestalHigh = cms.string(options.pedestalsHighGain)
process.hgcaltbrechits.gainLow = cms.string('')
process.hgcaltbrechits.gainHigh = cms.string('')

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")


process.output = cms.OutputModule("PoolOutputModule",
			fileName = cms.untracked.string(options.output)
                                 )

process.TFileService = cms.Service("TFileService", fileName = cms.string("Output.root"))

if(options.configuration == "1"):
    process.BadSpillFilter.layers_config = cms.int32(1)
    process.LayerSumAnalyzer.CERN_8layers_config = cms.int32(1)
    process.hgcaltbtrackingexample.CERN_8layers_config = cms.untracked.int32(0)
    process.hgcaltbshower.CERN_8layers_config = cms.untracked.int32(0)
elif(options.configuration == "2"):
    process.BadSpillFilter.layers_config = cms.int32(2)
    process.LayerSumAnalyzer.CERN_8layers_config = cms.int32(2)
    process.hgcaltbtrackingexample.CERN_8layers_config = cms.untracked.int32(1)
    process.hgcaltbshower.CERN_8layers_config = cms.untracked.int32(1)

#process.hgcaltbrechits.adcSaturation=cms.int32(1800)
#process.hgcaltbrechits.LG2HG=cms.double(9.75)

process.hgcaltbtrackingexample.minMip = cms.untracked.int32(9)
process.hgcaltbtrackingexample.maxMip = cms.untracked.int32(60)
process.hgcaltbtrackingexample.PrepareTreeForDisplay = cms.untracked.bool(False)
process.hgcaltbtrackingexample.doTrackCleaning = cms.untracked.bool(True)
process.hgcaltbtrackingexample.CMThreshold = cms.untracked.int32(100)

process.hgcaltbeventdisplay.minEnergy = cms.untracked.double(100)
process.hgcaltbeventdisplay.CMThreshold = cms.untracked.int32(100)

if (options.chainSequence == 1):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbdigisplotter)
elif (options.chainSequence == 2):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbtrackingexample)
elif (options.chainSequence == 3):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.hgcaltbeventdisplay)
elif (options.chainSequence == 4):
    process.hgcaltbcalotracks.maxEnergy=1e6
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.hgcaltbcalotracks*process.hgcaltbshower)
elif (options.chainSequence == 5):
    process.hgcaltbrechits.doCommonMode=False
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.myhgcaltbrechitsplotter)
elif (options.chainSequence == 6):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.LayerSumAnalyzer)
elif (options.chainSequence == 7):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.hgcaltbntuple)
elif (options.chainSequence == 8):
    process.hgcaltbcalotracks.doTrackCleaning=True
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.hgcaltbcalotracks*process.hgcaltbtrackanalyzer)

process.end = cms.EndPath(process.output)
