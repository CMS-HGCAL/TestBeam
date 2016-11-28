import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

repoFolder = "/afs/cern.ch/user/t/tquast/CMSSW_8_0_0_pre5/src/HGCal/"

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 #'/afs/cern.ch/user/t/tquast/eos/cms/store/group/upgrade/HGCAL/TestBeam/CERN/Sept2016/',       #use this for eos
                 '/user/data/Testbeam/September2016',        #use this for running on pclcd
                 #'/home/home1/institut_3a/quast/TestBeamSeptember2016/',        #use this for running on lx3a03
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw text input')

options.register('outputFolder',
                 '~/outputs/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Result of processing')

options.register('outputPostfix',
                 'default',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Postfix to the output file')
# options.register('commonPrefix',
#                  '',
#                  VarParsing.VarParsing.multiplicity.singleton,
#                  VarParsing.VarParsing.varType.string,
#                  'Input file to process')

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
                 repoFolder+'CondObjects/data/Ped_HighGain_L8.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to high gain pedestals file')

options.register('pedestalsLowGain',
                 repoFolder+'CondObjects/data/Ped_LowGain_L8.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to low gain pedestals file')

options.register('configuration',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 '-1 ADCtoMIP CERN; 0 ADCtoMIP FNAL; 1 if 8Layers with 5X0 sampling the center of the shower only; 2 if 8Layers with 25X0 sampling up to the tail of the shower')

options.register('reportEvery',
                100,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                'Frequency of event count print outs on the console')


'''Specific options for the position resolution'''
options.register('considerationMethod',
                 'all',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Possible arguments are: all, closest7, closest19, clusters'
                )
options.register('weightingMethod',
                 'squaredWeighting',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Possible arguments are: squaredWeighting, linearWeighting, logWeighting_5.0_1.0, logWeighting_5.0_0.5, logWeighting_7.0_1.0 '
                )
options.register('pedestalThreshold',
                 2.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Threshold for the pedestal subtraction calculation in the position resolution plugin. Unit is MIP.'
                )


options.maxEvents = -1

options.parseArguments()


# print options # <--- A better option is to call test_cfg is called with "print" argument e.g. cmsRun test_cfg.py print dataFolder=example1 outputFolder=example2 ...

if not os.path.isdir(options.dataFolder):
    print options.dataFolder
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
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

####################################
process.load('HGCal.StandardSequences.RawToDigi_cff')
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.dqm_cff')

process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(options.runNumber), ### maybe this should be read from the file
                            runEnergyMapFile = cms.untracked.string(repoFolder+"CondObjects/data/runEnergies.txt"), #the runs from the runEnergyMapFile are automatically added to the fileNames   
                            inputPathFormat=cms.untracked.string("file:%s/%s_Output_<RUN>.txt"%(options.dataFolder,options.runType)),  
                            fileNames=cms.untracked.vstring(["file:DUMMY"]), #'file:DUMMY'-->only files in the runEnergyMapFile are conidered
                                #["file:%s/%s_Output_%06d.txt"%(options.dataFolder,options.runType,options.runNumber) ]), ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            nSpills=cms.untracked.uint32(options.nSpills),
                            )



######
process.BadSpillFilter.configFile1 = "%s/CondObjects/data/Bad_Run_Spill_CFG1.txt" % repoFolder
process.BadSpillFilter.configFile2 = "%s/CondObjects/data/Bad_Run_Spill_CFG2.txt" % repoFolder
                           

process.hgcaltbdigisplotter.pedestalsHighGain = cms.untracked.string(options.pedestalsHighGain)
process.hgcaltbdigisplotter.pedestalsLowGain  = cms.untracked.string(options.pedestalsLowGain)

process.hgcaltbrechits.pedestalLow = cms.string(options.pedestalsLowGain)
process.hgcaltbrechits.pedestalHigh = cms.string(options.pedestalsHighGain)
process.hgcaltbrechits.gainLow = cms.string('')
process.hgcaltbrechits.gainHigh = cms.string('')

process.position_resolution_analyzer.considerationMethod = cms.string(options.considerationMethod)
process.position_resolution_analyzer.weightingMethod = cms.string(options.weightingMethod)
process.position_resolution_analyzer.pedestalThreshold = cms.double(options.pedestalThreshold)
process.position_resolution_analyzer.EventsFor2DGraphs = []    #first occuring events with that id are being documented with 2DGraphs 


process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")


if (options.chainSequence == 3):
    options.output = "%s/RECO_type%s_run%06d.root"%(options.outputFolder,options.runType,options.runNumber)
    process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string(options.output)
                                      )


# process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Display.root") )
if (options.chainSequence == 1):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Digi.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 3):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Unpacker_Check.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 4):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Reco_EventDisplay.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 5):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Reco.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 6):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Reco_Cluster.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 7):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Display_Cluster.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 8):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_Position_Resolution_%s.root"%(options.outputFolder,options.runType,options.outputPostfix)))


if(options.configuration == "-1"):
    process.BadSpillFilter.layers_config = cms.int32(-1)
    process.LayerSumAnalyzer.layers_config = cms.int32(-1)
    process.hgcaltbrechits.layers_config = cms.int32(-1)
elif(options.configuration == "0"):
    process.BadSpillFilter.layers_config = cms.int32(0)
    process.LayerSumAnalyzer.layers_config = cms.int32(0)
    process.hgcaltbrechits.layers_config = cms.int32(0)
elif(options.configuration == "1"):
    process.BadSpillFilter.layers_config = cms.int32(1)
    process.LayerSumAnalyzer.layers_config = cms.int32(1)
    process.hgcaltbrechits.layers_config = cms.int32(1)
    process.position_resolution_analyzer.layers_config = cms.int32(1)
elif(options.configuration == "2"):
    process.BadSpillFilter.layers_config = cms.int32(2)
    process.LayerSumAnalyzer.layers_config = cms.int32(2)
    process.hgcaltbrechits.layers_config = cms.int32(2)
    process.position_resolution_analyzer.layers_config = cms.int32(2)

########Activate this to produce event displays#########################################
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)

################Not needed for DQM purposes, produces digi histograms for each channel, and the pedestal txt file needed for Digi->Reco
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)

################Produces Reco histograms for each channel as well as a scatter plot of the Reco per channel#############
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)

#################Produces Clusters of Recos(7cells, 19cells and all cells(full hexagons only))################
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.LayerSumAnalyzer)

################Miscellaneous##############################################################################
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.FourLayerRecHitPlotterMax)

#Using chain sequence 3 only for testing purposes.
if (options.chainSequence == 1):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)
elif (options.chainSequence == 3):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits)
elif (options.chainSequence == 4):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)
elif (options.chainSequence == 5):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)
elif (options.chainSequence == 6):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.LayerSumAnalyzer)
elif (options.chainSequence == 7):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.hgcaltbeventdisplay)
elif (options.chainSequence == 8):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbclusters*process.position_resolution_analyzer)


# example for running display :
# cmsRun test_cfg_newEB.py runNumber=1291 runType=HGCRun nSpills=1 dataFolder='./' pedestalsHighGain="./CondObjects/data/pedHighGain1200.txt" pedestalsLowGain="./CondObjects/data/pedLowGain1200.txt" chainSequence=7 maxEvents=10


if (options.chainSequence == 3):
    process.end = cms.EndPath(process.output)

