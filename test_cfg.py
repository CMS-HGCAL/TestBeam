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
                 '0: if runType is PED then do Digi, if runType is HGC_Run then do Digi and Reco; 1: do Digi, 3: do both Digi and Reco (only Reco not implemented so far), 4: do event display sequence')

options.register('nSpills',
                 15,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of spills in run')

options.register('pedestalsHighGain',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to high gain pedestals file')

options.register('pedestalsLowGain',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to low gain pedestals file')

options.parseArguments()

if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

if not os.path.isdir(options.outputFolder):
    os.system("mkdir -p " + options.outputFolder)

if (options.runType != "PED" and options.runType != "HGCRun"):
    print options
    sys.exit("Error: only runtypes PED and HGCRun supported for now; given runType was %s"%(options.runType))

if (options.pedestalsHighGain == "" or options.pedestalsLowGain == ""): # Poor coding practice, but will have to do for now
    if (options.runType == "HGCRun"):
        options.pedestalsHighGain="CondObjects/data/Ped_HighGain_OneLayer_H2CERN.txt" #%(options.cmsswSource)
        options.pedestalsLowGain="CondObjects/data/Ped_LowGain_OneLayer_H2CERN.txt"  #%(options.cmsswSource)
    elif (options.runType == "PED"):
        options.pedestalsHighGain="%s/Ped_HighGain_OneLayer_%06d.txt"%(options.outputFolder, options.runNumber)
        options.pedestalsLowGain="%s/Ped_LowGain_OneLayer_%06d.txt"%(options.outputFolder, options.runNumber)

process = cms.Process("unpack")
process.load('HGCal.RawToDigi.hgcaltbdigis_cfi')
process.load('HGCal.RawToDigi.hgcaltbdigisplotter_cfi')
process.load('HGCal.Reco.hgcaltbrechitproducer_cfi')
process.load('HGCal.Reco.hgcaltbrechitplotter_cfi')

process.source = cms.Source("HGCalTBTextSource",
                            run=cms.untracked.int32(options.runNumber), ### maybe this should be read from the file
                            #fileNames=cms.untracked.vstring("file:Raw_data_New.txt") ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            fileNames=cms.untracked.vstring("file:%s/%s_Output_%06d.txt"%(options.dataFolder,options.runType,options.runNumber)), ### here a vector is provided, but in the .cc only the first one is used TO BE FIXED
                            nSpills=cms.untracked.uint32(options.nSpills)
)

process.hgcaltbdigisplotter = cms.EDAnalyzer("DigiPlotter",
                                     pedestalsHighGain=cms.untracked.string(options.pedestalsHighGain),
                                     pedestalsLowGain=cms.untracked.string(options.pedestalsLowGain)
                                     )

process.dumpRaw = cms.EDAnalyzer("DumpFEDRawDataProduct",
                              dumpPayload=cms.untracked.bool(True))

process.dumpDigi = cms.EDAnalyzer("HGCalDigiDump")


process.output = cms.OutputModule("PoolOutputModule",
			fileName = cms.untracked.string("test_output.root")
                                 )
process.hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                                OutputCollectionName = cms.string(''),
                                digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string(options.pedestalsLowGain),
                                pedestalHigh = cms.string(options.pedestalsHighGain),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                              )

# process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Display.root") )
if (options.chainSequence == 1):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Digi.root"%(options.outputFolder,options.runType,options.runNumber)))
elif (options.chainSequence == 3 or options.chainSequence == 4):
    process.TFileService = cms.Service("TFileService", fileName = cms.string("%s/%s_Output_%06d_Reco.root"%(options.outputFolder,options.runType,options.runNumber)))
# process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Layer.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("HGC_Output_6_Reco_Cluster.root") )


########Activate this to produce event displays#########################################
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)

################Not needed for DQM purposes, produces digi histograms for each channel, and the pedestal txt file needed for Digi->Reco
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)

################Produces Reco histograms for each channel as well as a scatter plot of the Reco per channel#############
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)

#################Produces Clusters of Recos(7cells, 19cells and all cells(full hexagons only))################
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.LayerSumAnalyzer)

if (options.chainSequence == 1):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbdigisplotter)
elif (options.chainSequence == 3):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)
    # process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter)
elif (options.chainSequence == 4):
    process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_new)
# process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.hgcaltbrechitsplotter_highgain_correlation_cm)
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.FourLayerRecHitPlotterMax)
#process.p =cms.Path(process.hgcaltbdigis*process.hgcaltbrechits*process.LayerSumAnalyzer)
process.end = cms.EndPath(process.output)
