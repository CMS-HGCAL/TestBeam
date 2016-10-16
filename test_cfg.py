import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'

options.register('dataFolder',
                 '/afs/cern.ch/work/r/rslu/public/HGC_TB_data_Sep2016/',
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
                 'CondObjects/data/Ped_HighGain_L8.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to high gain pedestals file')

options.register('pedestalsLowGain',
                 'CondObjects/data/Ped_LowGain_L8.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to low gain pedestals file')


options.register ('particle',
                  'electron',
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "particle type")


options.register ('energy',
                  '20GeV',
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "energy of beam")


## This output is disabled since we do not need edm root files (Eiko)
options.output = "test_output.root"

options.parseArguments()

# print options # <--- A better option is to call test_cfg is called with "print" argument e.g. cmsRun test_cfg.py print dataFolder=example1 outputFolder=example2 ...

if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

directory=options.particle+"/"+options.energy
if not os.path.isdir(directory):
    os.system("mkdir -p " + directory)

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


process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(options.runNumber), ### maybe this should be read from the file
                            #fileNames=cms.untracked.vstring("file:Raw_data_New.txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            fileNames=cms.untracked.vstring("file:%s/%s_Output_%06d.txt"%(options.dataFolder,options.runType,options.runNumber)), ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            nSpills=cms.untracked.uint32(options.nSpills),
)



######
process.hgcaltbdigisplotter.pedestalsHighGain = cms.untracked.string(options.pedestalsHighGain)
process.hgcaltbdigisplotter.pedestalsLowGain  = cms.untracked.string(options.pedestalsLowGain)

process.hgcaltbrechits.pedestalLow = cms.string(options.pedestalsLowGain)
process.hgcaltbrechits.pedestalHigh = cms.string(options.pedestalsHighGain)
process.hgcaltbrechits.gainLow = cms.string('')
process.hgcaltbrechits.gainHigh = cms.string('')

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")

process.LayerSumAnalyzer.particleType = cms.string(options.particle)


process.output = cms.OutputModule("PoolOutputModule",
			fileName = cms.untracked.string(options.output)
                                 )

if (options.chainSequence == 1):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Digi.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 3 or options.chainSequence == 4 or options.chainSequence == 5):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_%s_Output_%06d_Reco.root"%(directory,options.particle+options.energy,options.runType,options.runNumber)))


if (options.chainSequence == 1):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)
elif (options.chainSequence == 3):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)
elif (options.chainSequence == 4): 
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.LayerSumAnalyzer)
elif (options.chainSequence == 5):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm*process.hgcaltbrechitsplotter_highgain_new)


## Disabled since we do not need EDM root files (Eiko)
#process.end = cms.EndPath(process.output)
