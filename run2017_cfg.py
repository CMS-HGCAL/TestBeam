import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/disk2_2TB/July2017_TB_data_orm',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('runNumber',
                 106,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('outputFolder',
                 '/afs/cern.ch/work/a/asteen/public/data/july2017/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('evaluatePedestal',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Boolean to set for pedestal evaluation')

options.maxEvents = -1
options.output = "cmsswEvents.root"

options.parseArguments()
print options
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

electronicMap="HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt"
pedestalToCreateHighGain=options.outputFolder+"/pedestalHG_"+str(options.runNumber)+".txt"
pedestalToCreateLowGain=options.outputFolder+"/pedestalLG_"+str(options.runNumber)+".txt"
# pedestalToSubtractHighGain="pedestalHG_"+str(options.runNumber)+".txt"
# pedestalToSubtractLowGain="pedestalLG_"+str(options.runNumber)+".txt"
pedestalToSubtractHighGain="pedestalHG_1216.txt"
pedestalToSubtractLowGain="pedestalLG_1216.txt"


################################
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################

process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:%s/HexaData_Run%04d.raw"%(options.dataFolder,options.runNumber)),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30787),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(8),
                            NumberOfBytesForTheTrailer=cms.untracked.uint32(4),
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(4),
                            NumberOfOrmBoards=cms.untracked.uint32(1),
                            NSkipEvents=cms.untracked.uint32(0),
                            ReadTXTForTiming=cms.untracked.bool(False)
)

if options.evaluatePedestal==True:
    filename = options.outputFolder+"/PedestalOutput_"+str(options.runNumber)+".root"
    process.TFileService = cms.Service("TFileService", fileName=cms.string(filename))
else:
    filename = options.outputFolder+"/HexaOutput_"+str(options.runNumber)+".root"
    process.TFileService = cms.Service("TFileService", fileName = cms.string(filename))

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output)
)



process.content = cms.EDAnalyzer("EventContentAnalyzer") #add process.content in cms.Path if you want to check which collections are in the event

process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(True),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(electronicMap),
                                         HighGainPedestalFileName=cms.untracked.string(pedestalToCreateHighGain),
                                         LowGainPedestalFileName=cms.untracked.string(pedestalToCreateLowGain)
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
                                        HighGainPedestalFileName=cms.string(pedestalToSubtractHighGain),
                                        LowGainPedestalFileName=cms.string(pedestalToSubtractLowGain)
)

process.rawhitplotter = cms.EDAnalyzer("RawHitPlotter",
                                       InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       SubtractCommonMode=cms.untracked.bool(True)
)


if options.evaluatePedestal==True:
    process.p = cms.Path( process.rawdataplotter*process.pedestalplotter )
else:
    process.p = cms.Path( process.rawhitproducer*process.rawhitplotter )

#process.p = cms.Path( process.rawhitproducer*process.rawhitplotter*process.pulseshapeplotter )

process.end = cms.EndPath(process.output)

