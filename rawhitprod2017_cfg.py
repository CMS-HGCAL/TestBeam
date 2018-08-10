import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys,re

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('inputFileName',
                 './cmsswEvents.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('outputFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('electronicMap',
                 'HGCal/CondObjects/data/map_CERN_Hexaboard_OneModule_V2.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')
options.register('hgcalLayout',
                 'HGCal/CondObjects/data/layerGeom_1layer.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to hgcal layout file')

options.maxEvents = -1

options.parseArguments()
options.output=options.inputFileName
print options

if not os.path.isfile(options.inputFileName):
    sys.exit("Error: Input file not found or inaccessible!")



pedestalHighGain="pedestalHG.txt"
pedestalLowGain="pedestalLG.txt"
noisyChannels="noisyChannels.txt"

################################
process = cms.Process("rawhitprod")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s"%(options.inputFileName))
)

outputFileName=os.path.basename(options.inputFileName)
outputFileName,extent=os.path.splitext(outputFileName)
outputFileName=options.outputFolder+re.split('_edm',outputFileName)[0]+"_rawhit.root"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFileName)
)
process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_skiroc2cmsdata_*',
                                                                         'keep *_*_HGCALTBRAWHITS_*')
)

process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        ElectronicMap=cms.untracked.string(options.electronicMap),
                                        SubtractPedestal=cms.untracked.bool(True),
                                        MaskNoisyChannels=cms.untracked.bool(True),
                                        HighGainPedestalFileName=cms.untracked.string(pedestalHighGain),
                                        LowGainPedestalFileName=cms.untracked.string(pedestalLowGain),
                                        ChannelsToMaskFileName=cms.untracked.string(noisyChannels)
)

process.rawhitplotter = cms.EDAnalyzer("RawHitPlotter",
                                       InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                       ElectronicMap=cms.untracked.string(options.electronicMap),
                                       DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       SubtractCommonMode=cms.untracked.bool(True)
)

process.pulseshapeplotter = cms.EDAnalyzer("PulseShapePlotter",
                                           InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                           ElectronicMap=cms.untracked.string(options.electronicMap)
)


process.p = cms.Path( process.rawhitproducer*process.rawhitplotter*process.pulseshapeplotter )

process.end = cms.EndPath(process.output)
