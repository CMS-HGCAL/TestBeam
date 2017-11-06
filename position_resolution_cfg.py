import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys


options = VarParsing.VarParsing('standard')

####################################
# Options for reading in the data
options.register('repoFolder',
                '/afs/cern.ch/user/t/tquast/CMSSW_8_0_0_pre5/src/HGCal',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Directory to the repository to correctly navigate to some file in CondObjects'
                )

options.register('dataFolder',
                 '/home/data/Testbeam/September2016',        #use this for running on pclcd
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Main directory containing raw text input')

options.register('pathToRunEnergyFile',
                 '/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 '(Absolute) path to the file that indicates the runs to run the analysis on.')

options.register('readOnlyRuns',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.int,
                 'Run the analysis only for the indicated runs.'
                )

options.register('isData',
                 True,
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


options.register('alignmentParameterFiles',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Alignment parameter files as obtained from the (mille)pede framework.'
                )


####################################
#Specific options for the extra(inter)polation of impact positions in the MillepedeBinaryWrite and positionResolutionAnalyzer
options.register('useMWCReference',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Are the multi wire chamber information used for the exrapolation of referemce impact points?')

options.register('fittingMethod',
                'lineAnalytical',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Model of the electron tracks.'    
    )


####################################
#Specific options for reconstruction of impact positions in the MillepedeBinaryWrite and positionResolutionAnalyzer
options.register('considerationMethod',
                 'closest19',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Possible arguments are: all, closest7, closest19, clusters.'
                )

options.register('weightingMethod',
                 'logWeighting_3.5_1.0',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Possible arguments are: squaredWeighting, linearWeighting, logWeighting_<a>_1.0, .... '
                )

options.register('pedestalThreshold',
                 2.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Threshold for the common mode noise subtraction. Unit is MIP.'
                )

####################################
# Output file options
options.register('outputFolder',
                 '/home/outputs/Testbeam/September2016/',
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
if options.isData:
    process.source = cms.Source("HGCalTBTextSource",
                                runEnergyMapFile = cms.untracked.string(options.pathToRunEnergyFile), 
                                inputPathFormat=cms.untracked.string("file:%s/HGCRun_Output_<RUN>.txt"%options.dataFolder),  
                                MWCInputPathFormat=cms.untracked.string("file:%s/MWC/WC_H4Run<RUN>.txt"%options.dataFolder),
                                mwcRotation=cms.untracked.double(270.),
                                mwc2DeltaX=cms.untracked.double(0.),
                                mwc2DeltaY=cms.untracked.double(0.),
                                fileNames=cms.untracked.vstring(["file:DUMMY"]), #'file:DUMMY'-->only files in the runEnergyMapFile are considered
                                readOnlyRuns=cms.untracked.vint32(options.readOnlyRuns),
                                nSpills=cms.untracked.uint32(options.nSpills)
                                )
else:
    process.source = cms.Source("HGCalTBGenSimSource",
                            OutputCollectionName = cms.string(''), 
                            e_mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                            runEnergyMapFile = cms.untracked.string(options.pathToRunEnergyFile), 
                            inputPathFormat=cms.untracked.string("file:%s/<ENERGY>GeV/TBGenSim_<RUN>.root"%(options.dataFolder)),  
                            fileNames=cms.untracked.vstring(["file:DUMMY"]), #'file:DUMMY'-->only files in the runEnergyMapFile are considered,
                            energyNoise=cms.double(0.0),  
                            energyNoiseResolution=cms.double(0.0),
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
# position_resolution_analyzer options
process.position_resolution_analyzer.alignmentParameterFiles = cms.vstring(options.alignmentParameterFiles)
process.position_resolution_analyzer.fittingMethod = cms.string(options.fittingMethod)
process.position_resolution_analyzer.considerationMethod = cms.string(options.considerationMethod)
process.position_resolution_analyzer.weightingMethod = cms.string(options.weightingMethod)
process.position_resolution_analyzer.pedestalThreshold = cms.double(options.pedestalThreshold)
process.position_resolution_analyzer.fitPointWeightingMethod = cms.string("none")
process.position_resolution_analyzer.useMWCReference = cms.bool(options.useMWCReference)

####################################
# millepede_binarywriter options
process.millepede_binarywriter.fittingMethod = cms.string(options.fittingMethod)
process.millepede_binarywriter.considerationMethod = cms.string(options.considerationMethod)
process.millepede_binarywriter.weightingMethod = cms.string(options.weightingMethod)
process.millepede_binarywriter.pedestalThreshold = cms.double(options.pedestalThreshold)
process.millepede_binarywriter.fitPointWeightingMethod = cms.string("none")
process.millepede_binarywriter.totalEnergyThreshold = -1000.
process.millepede_binarywriter.useMWCReference = cms.bool(options.useMWCReference)
process.millepede_binarywriter.binaryFile = "%s/HGCRun_MillepedeBinary_%s.bin"%(options.outputFolder,options.outputPostfix)


####################################
# Necessary redefinition of sources if the input data is not read with the HGCalTBTextSource plugin
if not options.isData:
    process.hgcaltbclusters.rechitCollection = cms.InputTag("source","","unpack")
    process.position_resolution_analyzer.HGCALTBRECHITS = cms.InputTag("source","","unpack" )
    process.millepede_binarywriter.HGCALTBRECHITS = cms.InputTag("source","","unpack" )

                              
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


####################################
# Define the output file using TFileService. The MillepedeBinaryWriter defines its own, non-ROOT file as output
if (options.chainSequence == 8):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/HGCRun_Output_Position_Resolution_%s.root"%(options.outputFolder,options.outputPostfix)))



####################################
# Setup the sequences. With simulated data, 'process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits' do not run but are replaced by the HGCalTBGenSimSource plugin that converts
# the root tuples to HGCalRechits.
if options.isData:
    if (options.chainSequence == 8):
        process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.position_resolution_analyzer)
    elif (options.chainSequence == 9):
        process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.millepede_binarywriter)
        
else:
    if (options.chainSequence == 8):
        process.p =cms.Path(process.hgcaltbclusters*process.position_resolution_analyzer)
    elif (options.chainSequence == 9):
        process.p =cms.Path(process.hgcaltbclusters*process.millepede_binarywriter)






