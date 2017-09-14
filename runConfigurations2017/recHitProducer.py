import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:


options.register('dataFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/rawhits/RAWHITS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('processedFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/rechits/RECHITS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('outputFile',
                 '/home/tquast/tb2017/rechits/rechits_output_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('performPulseFit',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Perform the pulse fits for the rechit reconstruction. Has priority over averaging. ')

options.register('performAveraging',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Perform the averaging for the rechit reconstruction. ')

options.register('reportEvery',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.maxEvents = -1

options.parseArguments()
print options


electronicMap="HGCal/CondObjects/data/map_CERN_Hexaboard_July_6Layers.txt"

################################
process = cms.Process("rechitproducer")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s"%options.dataFile)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
)
process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.processedFile),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*',
                                                                         'keep *_*_RunData_*')
)

process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                        performPulseFit = cms.untracked.bool(bool(options.performPulseFit)),
                                        performAveraging = cms.untracked.bool(bool(options.performAveraging)),
                                        LG2HG = cms.untracked.vdouble(9.0,9.0,9.0,9.0,9.0,
                                                                      9.0,9.0,9.0,9.0,9.0),
                                        TOT2LG = cms.untracked.vdouble(9.0,9.0,9.0,9.0,9.0,
                                                                       9.0,9.0,9.0,9.0,9.0),
                                        HighGainADCSaturation = cms.untracked.vdouble(1800.,1800.,1800.,1800.,1800.,
                                                                       1800.,1800.,1800.,1800.,1800.),
                                        LowGainADCSaturation = cms.untracked.vdouble(1500.,1500.,1500.,1500.,1500.,
                                                                       1500.,1500.,1500.,1500.,1500.),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        TimeSample3ADCCut = cms.untracked.double(0.)
)

process.rechitplotter = cms.EDAnalyzer("RecHitPlotter",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       MipThreshold=cms.untracked.double(200),
                                       NoiseThreshold=cms.untracked.double(20)
)


process.p = cms.Path( process.rechitproducer*process.rechitplotter )

process.end = cms.EndPath(process.output)
