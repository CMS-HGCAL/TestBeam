import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/afs/cern.ch/user/r/rchatter/eos/cms/store/group/upgrade/HGCAL/TestBeam/CERN/Sept2016/',#modify path appropriately to where eos is mounted
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw text input')

options.register('outputFolder',
                 '/tmp/',#Choose the output directly where you wish the output root files to be written
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
process.load('HGCal.StandardSequences.dqm_cff')

#print "root://eoscms.cern.ch//eos/cms/%s/%s_Output_%06d.txt"%(options.dataFolder,options.runType,options.runNumber)

process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(options.runNumber), ### maybe this should be read from the file
                            #fileNames=cms.untracked.vstring("file:Raw_data_New.txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXE
                            fileNames=cms.untracked.vstring(["file:%s/%s_Output_%06d.txt"%(options.dataFolder,options.runType,options.runNumber), "file:%s/%s_Output_%06d.txt"%(options.dataFolder,options.runType,options.runNumber)]), ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            nSpills=cms.untracked.uint32(options.nSpills),
                            inputPathFormat = cms.untracked.string(""),
                            MWCInputPathFormat = cms.untracked.string(""),
                            mwcRotation = cms.untracked.double(0.),
                            mwc2DeltaX = cms.untracked.double(0.),
                            mwc2DeltaY = cms.untracked.double(0.),
                            readOnlyRuns = cms.untracked.vint32([]),
                            runEnergyMapFile = cms.untracked.string("runEnergyMapFile")
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


if (options.chainSequence == 7):
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
elif(options.configuration == "2"):
    process.BadSpillFilter.layers_config = cms.int32(2)
    process.LayerSumAnalyzer.layers_config = cms.int32(2)
    process.hgcaltbrechits.layers_config = cms.int32(2)
else:
    sys.exit("Error: Configuarion % is not supported in the position resolution analysis" % options.configuration)

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
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter)
elif (options.chainSequence == 4):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)
elif (options.chainSequence == 5):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm*process.hgcaltbrechitsplotter_highgain_new)
elif (options.chainSequence == 6):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits*process.LayerSumAnalyzer)
elif (options.chainSequence == 7):
    process.p =cms.Path(process.hgcaltbdigis*process.BadSpillFilter*process.hgcaltbrechits)


if (options.chainSequence == 7):
    process.end = cms.EndPath(process.output)