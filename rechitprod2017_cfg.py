import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:
options.register('dataFolder',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_june/HGCalTBRawHit',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing HGCalTBRawHit input')

options.register('runNumber',
                 24,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('edmOutputFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where edm output are stored')

options.register('rechitOutputFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('electronicMap',
                 'HGCal/CondObjects/data/map_CERN_Hexaboard_June_28Sensors_28EELayers_V1.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')

options.register('hgcalLayout',
                 'HGCal/CondObjects/data/layerGeom_june2018_h2_28layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to hgcal layout file')

options.register('adcCalibrationsFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/calibration/ADCCalibration_28modules_20-09-2018.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Root file containing ADC calibration parmaters')

options.maxEvents = -1
options.output = "cmsswEvents_RecHit.root"

options.parseArguments()
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

options.output = options.edmOutputFolder + "run%06d.root"%(options.runNumber)


print options

################################
process = cms.Process("rechitprod")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s/run%06d.root"%(options.dataFolder,options.runNumber)),
                            skipEvents=cms.untracked.uint32(0)
)

filename = options.rechitOutputFolder+"/HexaOutput_"+str(options.runNumber)+".root"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(filename)
)
process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.output),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*')
)

process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag('rawhitproducer','HGCALTBRAWHITS'),
                                        ElectronicsMap = cms.untracked.string(options.electronicMap),
                                        DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                        ADCCalibrationsFile=cms.untracked.string("file:%s"%(options.adcCalibrationsFile)),
                                        ExpectedMaxTimeSample=cms.untracked.int32(3),
                                        MaxADCCut=cms.untracked.double(20),
)

process.rechitplotter = cms.EDAnalyzer("RecHitPlotter",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       ElectronicsMap = cms.untracked.string(options.electronicMap),
                                       DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       MipThreshold=cms.untracked.double(5.0),
                                       NoiseThreshold=cms.untracked.double(0.5)
)

process.showeranalyzer = cms.EDAnalyzer("ShowerAnalyzer",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       NoiseThreshold=cms.untracked.double(0.5)
)

process.p = cms.Path( process.rechitproducer*process.rechitplotter )
#process.p = cms.Path( process.rechitproducer )

process.end = cms.EndPath(process.output)
