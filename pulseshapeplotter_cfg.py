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

options.register('outputFolder',
                 './',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('electronicMap',
                 'HGCal/CondObjects/data/map_CERN_Hexaboard_June_28Sensors_28EELayers_V0.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')

options.register('hgcalLayout',
                 'HGCal/CondObjects/data/layerGeom_june2018_h2_28layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to hgcal layout file')

options.maxEvents = -1

options.parseArguments()
if not os.path.isdir(options.dataFolder):
    sys.exit("Error: Data folder not found or inaccessible!")

print options

################################
process = cms.Process("pulseshapeplotter")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("file:%s/run%06d.root"%(options.dataFolder,options.runNumber)),
                            skipEvents=cms.untracked.uint32(0)
)

filename = options.outputFolder+"/HexaOutput_"+str(options.runNumber)+"_noCalib.root"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(filename)
)
process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string("./toto.root"),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*')
)

process.pulseshapeplotter = cms.EDAnalyzer("PulseShapePlotter",
                                           InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                           ElectronicMap=cms.untracked.string(options.electronicMap),
                                           DetectorLayout=cms.untracked.string(options.hgcalLayout),
                                           ExpectedMaxTimeSample=cms.untracked.int32(3),
                                           MaxADCCut=cms.untracked.double(20),
                                           SavePulseShapes=cms.untracked.bool(False)
)
process.p = cms.Path( process.pulseshapeplotter )
process.end = cms.EndPath(process.output)
