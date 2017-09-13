import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys


options = VarParsing.VarParsing('standard')

####################################
# Options for reading in the data
options.register('chainSequence',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Chain sequence to steer which process is run.'
                )

options.register('reportEvery',
                50000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )


options.register('outputFile',
                '/eos/user/t/tquast/outputs/Testbeam/July2017/default.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to the output file.'
                )

options.register('writeMinimal',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Write minimal in the DWC NTupelizer (1:yes, 0:no).'
                )

options.register('performAlignment',
                0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Perform alignment (1:yes, 0:no).'
                )

options.register('alignmentFiles',
                '/tmp/millepede.res',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Path to the alignment files.'
                )

options.register('inputFiles',
                [''],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the input files.'
                )

options.register('timingFiles',
                [''],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the timing files.'
                )

options.register('sumTriggerTimeStamps',
                [],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Sum the trigger time stamps in the timing file.'
                )

options.register('skipFirstNEvents',
                [],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Skip the first N events in the timing file.'
                )

options.register('triggerCountOffsets',
                [],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Indicate where the trigger count starts in the timing file.'
                )

options.register('runTypes',
                [''],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Run types (e.g. 1: 100 GeV pions).'
                )
       
#Alignment specific:
options.register('Layers',
                [],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 '0: E, 1: D, 2: A, 3: ext'
                )

options.register('coordinateString',
                'x',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Coordinate for translational alignment of DWCs (x/y).'
                )

options.register('outputMillepedeFile',
                '/tmp/millepede.bin',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to the millepede binary file.'
                )
    

options.parseArguments()

            

################################
# Setting an upper limit for the events to be processed, e.g. for debugging
options.maxEvents = -1
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################



####################################
# Initialize the data read-in plugins
process.source = cms.Source("HGCalTBWireChamberSource",
    OutputCollectionName = cms.string("WireChambers"), 
    fileNames = cms.untracked.vstring(["file:%s" % file for file in options.inputFiles]),
    timingFileNames = cms.vstring(["%s" % file for file in options.timingFiles]),
    sumTriggerTimes = cms.vint32([sumTrigger for sumTrigger in options.sumTriggerTimeStamps]),
    skipFirstNEvents = cms.vint32([skipFirstNEvents for skipFirstNEvents in options.skipFirstNEvents]),
    triggerCountOffsets = cms.vint32([triggerCountOffset for triggerCountOffset in options.triggerCountOffsets]),
    runType = cms.vstring([runType for runType in options.runTypes]),
    wc_resolution = cms.untracked.double(0.6),
    performAlignment = cms.untracked.bool(bool(options.performAlignment)),
    alignmentParamaterFiles = cms.vstring(options.alignmentFiles) 
)


####################################
#Millepede binary writer 
process.millepede_binarywriter.binaryFile = cms.string(options.outputMillepedeFile)
process.millepede_binarywriter.Layers = cms.vint32(options.Layers)
process.millepede_binarywriter.Coordinate = cms.string(options.coordinateString)
process.millepede_binarywriter.MWCQualityCut = cms.bool(False)
process.millepede_binarywriter.makeTree = cms.untracked.bool(True)
process.millepede_binarywriter.MWCHAMBERS = cms.InputTag("source","WireChambers","unpack")
process.millepede_binarywriter.RUNDATA = cms.InputTag("source","RunData","unpack")
process.millepede_binarywriter.fittingMethod = cms.string("lineAnalytical")
             

#DWC NTupelizer
process.dwc_ntupelizer.writeMinimal = cms.bool(bool(options.writeMinimal))


#Wire chamber producer
process.wirechamberproducer.OutputCollectionName = cms.string("DelayWireChambers") 
process.wirechamberproducer.RUNDATA = cms.InputTag("source","RunData","unpack")
process.wirechamberproducer.inputFile = cms.string("")


####################################
#add skip event exception which might occur for simulated samples because the last event is not properly passed forward
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))
if (options.chainSequence == 1):
    process.p = cms.Path(process.dwc_ntupelizer)
if (options.chainSequence == 2):
    process.p = cms.Path(process.millepede_binarywriter)


#process.p = cms.Path(process.wirechamberproducer)
