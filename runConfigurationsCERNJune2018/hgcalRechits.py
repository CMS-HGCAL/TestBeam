import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/home/tquast/tbJune2018_H2/hgcalReco/rawhits_000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'file containing raw input')

options.register('processedFile',
                 '/home/tquast/tbJune2018_H2/hgcalReco/rechits_000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('outputFile',
                 '/home/tquast/tbJune2018_H2/hgcalReco/rechit_analysis_000223.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('NHexaBoards',
                28,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )


options.register('electronicMap',
                 "map_CERN_Hexaboard_June_28Sensors_28EELayers_V0.txt",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('hgcalLayout',
                 'layerGeom_june2018_h2_28layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')

options.register('adcCalibrations',
                 'hgcal_calibration_June2018_v0.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal ADC to MIP calibration file in HGCal/CondObjects/data/')

options.register('reportEvery',
                100,
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
                                        RUNDATA = cms.InputTag("source", "RunData" ), 
                                        GlobalTimestampCollectionName=cms.InputTag("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        MaskNoisyChannels=cms.untracked.bool(bool(0)),
                                        ChannelsToMaskFileName=cms.untracked.string(""),
                                        NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                        TimeSample3ADCCut = cms.untracked.double(15.),
                                        investigatePulseShape = cms.untracked.bool(True),
                                        timingNetworks = cms.untracked.string("")
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

process.p = cms.Path(process.rechitproducer)

process.end = cms.EndPath(process.output)
