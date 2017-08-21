import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/july2017/July2017_TB_data_orm',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 106,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('chainSequence',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Chain sequence to run.')

options.register('outputFile',
                 '/afs/cern.ch/work/a/asteen/public/data/july2017/HexaOutput_106.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('DWCFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/reconstructed_DWC_data_full.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the reconstructed DWC file.')

options.register('reportEvery',
                100,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                'Frequency of event count printouts on the console.')

options.register('evaluatePedestal',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Boolean to set for pedestal evaluation')

options.register('timingFile',
                '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the timing file which is produced with chain sequence 0.')

options.maxEvents = -1


options.output = "cmsswEvents.root"

options.parseArguments()
print options

#for the quick check for synchronisation
if options.chainSequence==2:
    options.maxEvents = -1
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

electronicMap="HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt"


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

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:%s/HexaData_Run%04d.raw"%(options.dataFolder,options.runNumber)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30787),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(8 if options.runNumber<1241 else 12),          #for the new headers/trailers from run 1241 onward: 12
                            NumberOfBytesForTheTrailer=cms.untracked.uint32(4),         #for the new headers/trailers from run 1241 onward: 4
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(4 if options.runNumber<1241 else 12),   #for the new headers/trailers from run 1241 onward: 12
                            NSkipEvents=cms.untracked.uint32(0),
                            ReadTXTForTiming=cms.untracked.bool(False),
                            timingFilePath=cms.untracked.string(options.timingFile),
                            timingFiles=cms.vstring("%s/HexaData_Run%04d_TIMING_RDOUT_ORM0.txt"%(options.dataFolder,options.runNumber),
                                                    "%s/HexaData_Run%04d_TIMING_RDOUT_ORM1.txt"%(options.dataFolder,options.runNumber),
                                                    "%s/HexaData_Run%04d_TIMING_RDOUT_ORM2.txt"%(options.dataFolder,options.runNumber))
                            )


process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output)
)

outputFilePath = "/".join(options.outputFile.split("/")[0:-1])

pedestalHighGain="%s/pedestalHG_%s.txt" % (outputFilePath, str(options.runNumber))
pedestalLowGain="%s/pedestalLG_%s.txt" % (outputFilePath, str(options.runNumber))

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))

process.content = cms.EDAnalyzer("EventContentAnalyzer") #add process.content in cms.Path if you want to check which collections are in the event


process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(True),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(electronicMap),
                                         HighGainPedestalFileName=cms.untracked.string(pedestalHighGain),
                                         LowGainPedestalFileName=cms.untracked.string(pedestalLowGain)
                                         )


process.rawdataplotter = cms.EDAnalyzer("RawDataPlotter",
                                        SensorSize=cms.untracked.int32(128),
                                        EventPlotter=cms.untracked.bool(False),
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata")
                                        )

process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        ElectronicMap=cms.untracked.string(electronicMap),
                                        SubtractPedestal=cms.untracked.bool(True),
                                        HighGainPedestalFileName=cms.string(pedestalHighGain),
                                        LowGainPedestalFileName=cms.string(pedestalLowGain)
                                        )

process.rawhitplotter = cms.EDAnalyzer("RawHitPlotter",
                                       InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                       MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers","unpack"),
                                       RUNDATA = cms.InputTag("source","RunData","unpack"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       SubtractCommonMode=cms.untracked.bool(True),
                                       outputFile = cms.untracked.string(options.outputFile)
                                       )

process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag('rawhitproducer','HGCALTBRAWHITS'),
                                        LG2HG = cms.untracked.vdouble(10.0),
                                        TOT2LG = cms.untracked.vdouble(10.0),
                                        HighGainADCSaturation = cms.untracked.double(2000),
                                        LowGainADCSaturation = cms.untracked.double(1500),
                                        ElectronicsMap = cms.untracked.string('HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt'),
                                        CommonModeThreshold = cms.untracked.double(100.),
                                        KeepOnlyTimeSample3 = cms.untracked.bool(False)
                                        )

#Wire chamber producer
process.wirechamberproducer.OutputCollectionName = cms.string("DelayWireChambers") 
process.wirechamberproducer.RUNDATA = cms.InputTag("source","RunData","unpack")
process.wirechamberproducer.inputFile = cms.string(options.DWCFile)


####################################
#add skip event exception which might occur for simulated samples because the last event is not properly passed forward
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


if options.chainSequence==0:
    process.p = cms.Path()
if options.chainSequence==1:
    process.p = cms.Path( process.rawdataplotter*process.pedestalplotter )
if options.chainSequence==2:
    process.p = cms.Path( process.rawhitproducer*process.wirechamberproducer*process.rawhitplotter)



#process.end = cms.EndPath(process.output)

