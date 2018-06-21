import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('processedFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('outputFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('NHexaBoards',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('electronicMap',
                 'map_DESY_March2018_config4_V0.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('adcCalibrations',
                 'hgcal_calibration.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal ADC to MIP calibration file in HGCal/CondObjects/data/')

options.register('hgcalLayout',
                 'layerGeom_desymarch2018_configuration4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')


options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 ''
                )

options.maxEvents = -1

options.parseArguments()
print options


electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap
hgcalLayout="HGCal/CondObjects/data/%s" % options.hgcalLayout
adcCalibrations="HGCal/CondObjects/data/%s" % options.adcCalibrations

################################
process = cms.Process("rechits")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################



process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))


process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s"%options.dataFile)
)

process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag("rawhitproducer", "HGCALTBRAWHITS"),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        ExpectedMaxTimeSample=cms.untracked.int32(3),
                                        MaxADCCut=cms.untracked.double(20)
)

process.rechitplotter = cms.EDAnalyzer("RecHitPlotter",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       MipThreshold=cms.untracked.double(200),
                                       NoiseThreshold=cms.untracked.double(20)
)



####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.processedFile),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*',
                                                                         'keep *_*_HGCALGLOBALTIMESTAMPS_*',
                                                                         'keep *_*_RunData_*')
)


process.p = cms.Path(process.rechitproducer*process.rechitplotter)

process.end = cms.EndPath(process.output)
