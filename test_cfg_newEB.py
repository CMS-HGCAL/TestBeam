import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'

options.register('dataFolder',
                 '~/eos/cms/store/group/upgrade/HGCAL/TestBeam/CERN/Sept2016/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw text input')

options.register('outputFolder',
                 '/tmp/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Result of processing')

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
                 'CondObjects/data/Ped_HighGain_L8.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to high gain pedestals file')

options.register('pedestalsLowGain',
                 'CondObjects/data/Ped_LowGain_L8.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to low gain pedestals file')

options.register('configuration',
                 -1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 '-1 ADCtoMIP CERN; 0 ADCtoMIP FNAL; 1 if 8Layers with 5X0 sampling the center of the shower only; 2 if 8Layers with 25X0 sampling up to the tail of the shower')



options.output = "test_output.root"
options.maxEvents = -1

options.parseArguments()

# print options # <--- A better option is to call test_cfg is called with "print" argument e.g. cmsRun test_cfg.py print dataFolder=example1 outputFolder=example2 ...

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
process.load('HGCal.StandardSequences.TrackingReco_cff')
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


process.output = cms.OutputModule("PoolOutputModule",
			fileName = cms.untracked.string(options.output)
                                 )

# process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Display.root") )
if (options.chainSequence == 1):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Digi.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 3):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Unpacker_Digi_Check.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 4):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Reco_EventDisplay.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 5):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Reco.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 6):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Reco_Cluster.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 7):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Display_Cluster.root"%(options.outputFolder,options.runType,options.runNumber)))



if(options.configuration == "-1"):
    process.LayerSumAnalyzer.layers_config = cms.int32(-1)
    process.hgcaltbrechits.layers_config = cms.int32(-1)
elif(options.configuration == "0"):
    process.LayerSumAnalyzer.layers_config = cms.int32(0)
    process.hgcaltbrechits.layers_config = cms.int32(0)
elif(options.configuration == "1"):
    process.LayerSumAnalyzer.layers_config = cms.int32(1)
    process.hgcaltbrechits.layers_config = cms.int32(1)
elif(options.configuration == "2"):
    process.LayerSumAnalyzer.layers_config = cms.int32(2)
    process.hgcaltbrechits.layers_config = cms.int32(2)

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

if (options.chainSequence == 1):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)
elif (options.chainSequence == 3):
    process.p =cms.Path(process.hgcaltbdigis)
elif (options.chainSequence == 4):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)
elif (options.chainSequence == 5):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)
elif (options.chainSequence == 6):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.LayerSumAnalyzer)
elif (options.chainSequence == 7):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbclusters*process.hgcaltbeventdisplay)
elif (options.chainSequence == 8):
    process.hgcaltbcalotracks.doTrackCleaning=True
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbclusters*process.hgcaltbcalotracks*process.hgcaltbtrackanalyzer)
# example for running display :
# cmsRun test_cfg_newEB.py runNumber=1291 runType=HGCRun nSpills=1 dataFolder='./' pedestalsHighGain="./CondObjects/data/pedHighGain1200.txt" pedestalsLowGain="./CondObjects/data/pedLowGain1200.txt" chainSequence=7 maxEvents=10
# example for running tracking :
# cmsRun test_cfg_newEB.py runNumber=1312 runType=HGCRun nSpills=15 dataFolder='./' pedestalsHighGain="./CondObjects/data/pedHighGain1200.txt" pedestalsLowGain="./CondObjects/data/pedLowGain1200.txt" chainSequence=8

process.end = cms.EndPath(process.output)
