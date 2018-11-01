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
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.register('outputFile',
                '~/tmp/test_512.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to the output file.'
                )

options.register('performAlignment',
                1,
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
                ['/home/tquast/tbOctober2018_H2/reco_VME/unpacked_VME_run000512.root'],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the input files.'
                )

options.register('timingFiles',
                ['/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/ORM_timingFiles/timing_512.txt'],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the timing files.'
                )

options.register('triggerTimeDifferenceTolerance',
                0.5,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Accepted tolerance in trigger time stamp difference for the synchronisation.'
                )

options.register('TDCTriggerTimeStampConversionToMs',
                8./10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Conversion from TDC trigger time stamp to ms'
                )

options.register('skipFirstNEvents',
                [1],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Skip the first N events in the timing file.'
                )

options.register('triggerCountOffsets',
                [582],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Indicate where the trigger count starts in the timing file.'
                )

options.register('skipTDCTriggers',
                [0],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Can events in the TDC be skipped?.'
                )

options.register('setupIDs',
                [22],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Setup IDs.'
                )

options.register('pdgIDs',
                [211],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'PDG IDs.'
                )

options.register('beamEnergies',
                [300.],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.float,
                 'beamEnergies [GeV].'
                )

options.register('areaSpecification',
                "H2_October2018",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Run types (e.g. 1: 100 GeV pions).'
                )

#Alignment specific:
options.register('Layers',
                [0],
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
process = cms.Process("dwcReco")
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
    sumTriggerTimes = cms.vint32([1]),
    skipFirstNEvents = cms.vint32([skipFirstNEvents for skipFirstNEvents in options.skipFirstNEvents]),
    triggerCountOffsets = cms.vint32([triggerCountOffset for triggerCountOffset in options.triggerCountOffsets]),
    triggerTimeDifferenceTolerance = cms.untracked.double(options.triggerTimeDifferenceTolerance),
    TDCTriggerTimeStampConversionToMs = cms.untracked.double(options.TDCTriggerTimeStampConversionToMs),
    allowForTDCEventSkipping = cms.vint32([skipTDCTrigger for skipTDCTrigger in options.skipTDCTriggers]),
    setupIDs = cms.vint32([setupID for setupID in options.setupIDs]),
    pdgIDs = cms.vint32([pdgID for pdgID in options.pdgIDs]),
    beamEnergies = cms.vdouble([energies for energies in options.beamEnergies]),
    triggerTimingFormat = cms.vint32([1]),
    hitsPerChannelStored = cms.vint32([1]*32),
    areaSpecification = cms.untracked.string(options.areaSpecification),
    wc_resolutions = cms.untracked.vdouble([1.0, 1.0, 1.0, 1.0]),
    performAlignment = cms.untracked.bool(bool(options.performAlignment)),
    alignmentParamaterFiles = cms.vstring(options.alignmentFiles) 
)

#DWC NTupelizer
if (options.chainSequence == 1):
    process.dwc_ntupelizer.MWCHAMBERS = cms.InputTag("source","WireChambers","dwcReco" )
    process.dwc_ntupelizer.RUNDATA = cms.InputTag("source","RunData","dwcReco" )
    process.dwc_ntupelizer.writeMinimal = cms.bool(False)


####################################
#Millepede binary writer 
if (options.chainSequence == 2):
    process.millepede_binarywriter.binaryFile = cms.string(options.outputMillepedeFile)
    process.millepede_binarywriter.Layers = cms.vint32(options.Layers)
    process.millepede_binarywriter.Coordinate = cms.string(options.coordinateString)
    process.millepede_binarywriter.MWCQualityCut = cms.bool(False)
    process.millepede_binarywriter.makeTree = cms.untracked.bool(True)
    process.millepede_binarywriter.MWCHAMBERS = cms.InputTag("source","WireChambers","dwcReco")
    process.millepede_binarywriter.RUNDATA = cms.InputTag("source","RunData","dwcReco")
    process.millepede_binarywriter.fittingMethod = cms.string("lineAnalytical")

            
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


