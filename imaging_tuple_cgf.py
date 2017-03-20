import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os, sys


options = VarParsing.VarParsing('standard')

####################################
# Options for reading in the data
options.register('repoFolder',
                '/afs/cern.ch/user/t/tquast/CMSSW_8_0_0_pre5/src/HGCal',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Directory to the repository to correctly navigate to some file in CondObjects'
                )

options.register('fileNames',
                 #'file:/home/data/MC/September2016/70GeV/TBGenSim_1.root',       
                 'file:/home/data/MC/EnergyReco/February2017/11/test.root',       
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Input files.')

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
                 10,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Here: only 10 is configured.')

options.register('reportEvery',
                100,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                'Frequency of event count printouts on the console.')


####################################
# Options related to the experimental setup
options.register('configuration',
                 "1",
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
if options.isData:
    process.source = cms.Source("HGCalTBTextSource",
                                runEnergyMapFile = cms.untracked.string(""), 
                                inputPathFormat=cms.untracked.string(""),  
                                MWCInputPathFormat=cms.untracked.string(""),
                                mwcRotation=cms.untracked.double(0.),
                                mwc2DeltaX=cms.untracked.double(0.),
                                mwc2DeltaY=cms.untracked.double(0.),
                                fileNames=cms.untracked.vstring(options.fileNames), 
                                readOnlyRuns=cms.untracked.vint32([]),
                                nSpills=cms.untracked.uint32(options.nSpills)
                                )
else:
    process.source = cms.Source("HGCalTBGenSimSource",
                            OutputCollectionName = cms.string(''), 
                            e_mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                            runEnergyMapFile = cms.untracked.string(""), 
                            inputPathFormat=cms.untracked.string(""),  
                            fileNames=cms.untracked.vstring(options.fileNames), 
                            energyNoise=cms.double(0.0),  
                            energyNoiseResolution=cms.double(0.0),
                            createMWC=cms.bool(False),
                            MWCSmearingResolution=cms.double(50.)     #value is in microns! 
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
# Necessary redefinition of sources if the input data is not read with the HGCalTBTextSource plugin
if not options.isData:
    process.imaging_tuple_writer.HGCALTBRECHITS = cms.InputTag("source","","unpack" )
    process.LayerSumAnalyzer.HGCALTBRECHITS = cms.InputTag("source","","unpack" )

                              
####################################
# Set the configuration in the appropriate plugins that need this information
if(options.configuration == "1"):
    process.BadSpillFilter.layers_config = cms.int32(1)
    process.hgcaltbrechits.layers_config = cms.int32(1)
    process.LayerSumAnalyzer.layers_config = cms.int32(1)
elif(options.configuration == "2"):
    process.BadSpillFilter.layers_config = cms.int32(2)
    process.hgcaltbrechits.layers_config = cms.int32(2)
    process.LayerSumAnalyzer.layers_config = cms.int32(2)
else:
    sys.exit("Error: Configuarion %s is not supported in the position resolution analysis" % options.configuration)


####################################
# Define the output file using TFileService. The MillepedeBinaryWriter defines its own, non-ROOT file as output
if (options.chainSequence == 10):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/HGCRun_Output_Imaging_Tuple_%s.root"%(options.outputFolder,options.outputPostfix)))
elif (options.chainSequence == 6):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/HGCRun_LayerSumAnalyzer_%s.root"%(options.outputFolder,options.outputPostfix)))


####################################
# Setup the sequences. With simulated data, 'process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits' do not run but are replaced by the HGCalTBGenSimSource plugin that converts
# the root tuples to HGCalRechits.
if options.isData:
    if (options.chainSequence == 10):
        process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.imaging_tuple_writer)
    elif (options.chainSequence == 6):
        process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.LayerSumAnalyzer)    

else:
    if (options.chainSequence == 10):
        process.p =cms.Path(process.imaging_tuple_writer)
    elif (options.chainSequence == 6):
        process.p =cms.Path(process.LayerSumAnalyzer)       






