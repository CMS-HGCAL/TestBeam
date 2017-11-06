import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys


options = VarParsing.VarParsing('standard')

####################################
# Options for reading in the data
options.register('repoFolder',
                '/afs/cern.ch/work/r/rchatter/CERN_TestBeam_2017/CMSSW_8_0_21/src/HGCal',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Directory to the repository to correctly navigate to some file in CondObjects'
                )

options.register('dataFolder',
                 '/eos/cms/store/group/upgrade/HGCAL/simulation/8moduleIv8_PhysicsList_FTFP_BERT_EMM_allImprovements_2017Jan/mc/CRAB_PrivateMC/crab_Ele100GeV/170102_203458/0000',        #use this for running on pclcd
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Main directory containing raw text input')

options.register('readOnlyRuns',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Run the analysis only for the indicated runs.'
                )

options.register('isData',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Is the analysis run on real data (otherwise on simulated samples)?')

options.register('nSpills',
                 15,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of spills per run before read-out is stopped.')


####################################
# Options specific for the pedestal subtraction
options.register('pedestalsHighGain',
                 '%s/CondObjects/data/Ped_HighGain_L8.txt' % options.repoFolder,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to high gain pedestals file.')

options.register('pedestalsLowGain',
                 '%s/CondObjects/data/Ped_LowGain_L8.txt' % options.repoFolder,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to low gain pedestals file.')

####################################
# Options for running the analysis, e.g. indicating what sequence is to be run
options.register('chainSequence',
                 8,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Here: only 8 (position resolution) and 9 (writing of millepede binary) is configured.')

options.register('reportEvery',
                100,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                'Frequency of event count printouts on the console.')

####################################
# Options related to the experimental setup
options.register('configuration',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 '1 if 8Layers with 5X0 sampling the center of the shower only; 2 if 8Layers with 25X0 sampling up to the tail of the shower')

####################################
# Output file options
options.register('outputFolder',
                 '/tmp/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Main directory of the output files.')

options.register('outputPostfix',
                 'default',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Postfix to the output file.')


options.parseArguments()
####################################
# Check if necessary folders exist
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

if not os.path.isdir(options.outputFolder):
    os.system("mkdir -p " + options.outputFolder)

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
process.load('HGCal.StandardSequences.RawToDigi_cff')
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.dqm_cff')
####################################


####################################
# Initialize the data read-in plugins
if not options.isData:
    process.source = cms.Source("HGCalTBGenSimSource_Only",
                            OutputCollectionName = cms.string(''),
                            e_mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                            fileNames=cms.untracked.vstring(["file:/eos/cms/store/group/upgrade/HGCAL/simulation/8moduleIv8_PhysicsList_FTFP_BERT_EMM_allImprovements_2017Jan/mc/CRAB_PrivateMC/crab_Ele100GeV/170102_203458/0000/TBGenSim_803.root"]) 
                            )


####################################
#add skip event exception which might occur for simulated samples because the last event is not properly passed forward
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

####################################
# BadSpillFilter options
process.BadSpillFilter.nameCFG1 = cms.string("%s/CondObjects/data/Bad_Run_Spill_CFG1.txt" % options.repoFolder)
process.BadSpillFilter.nameCFG2 = cms.string("%s/CondObjects/data/Bad_Run_Spill_CFG2.txt" % options.repoFolder)
                          
####################################
# hgcaltbdigisplotter options
process.hgcaltbdigisplotter.pedestalsHighGain = cms.untracked.string(options.pedestalsHighGain)
process.hgcaltbdigisplotter.pedestalsLowGain  = cms.untracked.string(options.pedestalsLowGain)

####################################
# hgcaltbrechits options
process.hgcaltbrechits.pedestalLow = cms.string(options.pedestalsLowGain)
process.hgcaltbrechits.pedestalHigh = cms.string(options.pedestalsHighGain)
process.hgcaltbrechits.gainLow = cms.string('')
process.hgcaltbrechits.gainHigh = cms.string('')

####################################
# Set the configuration in the appropriate plugins that need this information
if(options.configuration == "1"):
    process.BadSpillFilter.layers_config = cms.int32(1)
    process.hgcaltbrechits.layers_config = cms.int32(1)
    process.position_resolution_analyzer.layers_config = cms.int32(1)
elif(options.configuration == "2"):
    process.BadSpillFilter.layers_config = cms.int32(2)
    process.hgcaltbrechits.layers_config = cms.int32(2)
    process.position_resolution_analyzer.layers_config = cms.int32(2)
else:
    sys.exit("Error: Configuarion % is not supported in the position resolution analysis" % options.configuration)
######################################
process.output = cms.OutputModule("PoolOutputModule",
                        fileName = cms.untracked.string('Test_SimToTBRecHits_803.root')
                                 )
####################################
# Define the output file using TFileService. The MillepedeBinaryWriter defines its own, non-ROOT file as output
if (options.chainSequence == 8):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/HGCRun_Output_Position_Resolution_%s.root"%(options.outputFolder,options.outputPostfix)))

####################################
# Setup the sequences. With simulated data, 'process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits' do not run but are replaced by the HGCalTBGenSimSource plugin that converts
# the root tuples to HGCalRechits.
if not options.isData:
    if (options.chainSequence == 8):
        process.p =cms.Path()

process.end = cms.EndPath(process.output)
