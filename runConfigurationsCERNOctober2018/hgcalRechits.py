import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/rawhits/v1/RAWHITS_419.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'file containing raw input')

options.register('outputFile',
                 'toto.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('processedFile',
                 'rechit.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('NHexaBoards',
                94,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('electronicMap',
                 'emap_full_October2018_setup1_v5_promptReco.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('hgcalLayout',
                 'layer_geom_full_October2018_setup1_v5_promptReco.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')

options.register('adcCalibrations',
                 'hgcal_calibration_October2018_JuneCalibForEE.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal ADC to MIP calibration file in HGCal/CondObjects/data/')

options.register('ExpectedMaxTimesample',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of the timesample where the maximum of the pulse is expected.'
                )

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
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        calibrationPerChannel=cms.untracked.bool(True),
                                        ExpectedMaxTimeSample=cms.untracked.int32(options.ExpectedMaxTimesample),
                                        MaxADCCut=cms.untracked.double(15),
                                        subtractCommonMode=cms.untracked.bool(True),
                                        subtractCommonModeOption=cms.untracked.string("board_thr"),
                                        commonModeThreshold=cms.untracked.double(100.),
                                        TSForCommonModeNoiseSubtraction=cms.untracked.int32(-1)
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
