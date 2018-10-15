import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'

options.register('runNumber',
                 313,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Input run to process')

options.register('beamEnergy',
                30,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Beam energy.'
                )

options.register('beamParticlePDGID',
                11,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Beam particles PDG ID.'
                )






# no need to change anything from here below 

options.register('dataFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/unpacked/run%06d.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('timingFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/ORM_timingFiles/timing_%s.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('outputFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/prompt_reco/v2/prompt_reco_%s.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')


options.register('runType',
                 "Beam",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Run type: Pedestal, Beam, Simulation.'
                )

options.register('setupConfiguration',
                22,     #22: October 2018
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'setupConfiguration (1: July - 4: 20 Layers in October in H6A".'
                )

options.register('pedestalHighGainFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/pedestalFiles/pedestalHG_259.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('pedestalLowGainFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/pedestalFiles/pedestalLG_259.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('noisyChannelsFile',
                 '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/pedestalFiles/noisyChannels_259.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('electronicMap',
                 "emap_full_October2018_v3_promptReco.txt",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'path to the electronic map')

options.register('NHexaBoards',
                94,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('NLayers',
                40,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of layers for analysis.'
                )

options.register('hgcalLayout',
                 'layer_geom_full_October2018_v3_promptReco.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')

options.register('adcCalibrations',
                 'hgcal_calibration_October2018_v0_promptReco.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal ADC to MIP calibration file in HGCal/CondObjects/data/')

options.register('layerPositionFile',
                 '/afs/cern.ch/user/d/daq/CMSSW_9_3_0/src/HGCal/CondObjects/data/layer_distances_CERN_Hexaboard_October2018_promptReco_v0.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('SubtractPedestal',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Subtract the pedestals.'
                )

options.register('MaskNoisyChannels',
                0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.register('reportEvery',
                10,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 '.'
                )

options.maxEvents = -1

options.parseArguments()
print options

electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap
hgcalLayout="HGCal/CondObjects/data/%s" % options.hgcalLayout
adcCalibrations="HGCal/CondObjects/data/%s" % options.adcCalibrations
layerPositionFile=options.layerPositionFile

################################
process = cms.Process("CMSSWUnpacker")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

################################
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile % options.runNumber))

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')


process.source = cms.Source("HGCalTBEUDAQDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:%s"%options.dataFile % options.runNumber),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NSkipEvents=cms.untracked.uint32(0),
                            runNumber=cms.untracked.int32(options.runNumber),
                            beamEnergy=cms.untracked.double(options.beamEnergy),
                            beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),
                            runType=cms.untracked.string(options.runType),
                            setupConfiguration=cms.untracked.uint32(options.setupConfiguration)
)

process.timingfilewriter = cms.EDAnalyzer("HGCalTBTimingFileWriter",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        TimingFilePath=cms.untracked.string(options.timingFile % options.runNumber)
)

process.pedestalplotter = cms.EDAnalyzer("PedestalPlotter",
                                         SensorSize=cms.untracked.int32(128),
                                         WritePedestalFile=cms.untracked.bool(False),
                                         InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                         ElectronicMap=cms.untracked.string(electronicMap),
                                         HighGainPedestalFileName=cms.untracked.string(options.pedestalHighGainFile),
                                         LowGainPedestalFileName=cms.untracked.string(options.pedestalLowGainFile),
                                         WriteNoisyChannelsFile=cms.untracked.bool(False),
                                         NoisyChannelsFileName=cms.untracked.string(options.noisyChannelsFile),
                                         NTSForPedestalComputation=cms.untracked.int32(0)
)

process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        GlobalTimestampCollectionName=cms.string("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicMap=cms.untracked.string(electronicMap),
                                        SubtractPedestal=cms.untracked.bool(bool(options.SubtractPedestal)),
                                        MaskNoisyChannels=cms.untracked.bool(bool(options.MaskNoisyChannels)),
                                        HighGainPedestalFileName=cms.untracked.string(options.pedestalHighGainFile),
                                        LowGainPedestalFileName=cms.untracked.string(options.pedestalLowGainFile),
                                        ChannelsToMaskFileName=cms.untracked.string(options.noisyChannelsFile)
)


process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag("rawhitproducer", "HGCALTBRAWHITS"),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        calibrationPerChannel=cms.untracked.bool(True),
                                        ExpectedMaxTimeSample=cms.untracked.int32(3),
                                        MaxADCCut=cms.untracked.double(15)
)


process.eventdisplay = cms.EDAnalyzer("EventDisplay",
                                RUNDATA = cms.InputTag("source", "RunData"), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                electronicsMap = cms.untracked.string(electronicMap),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                eventsToPlot=cms.vint32(range(1, 4))
                              )


process.rechitntupler = cms.EDAnalyzer("RecHitNtupler",
                                       InputCollection=cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                       RUNDATA = cms.InputTag("source", "RunData"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       layerPositionFile = cms.untracked.string(layerPositionFile),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(True),
                                       MipThreshold=cms.untracked.double(2.0),
                                       NoiseThreshold=cms.untracked.double(0.0)
)


process.variablecomputation = cms.EDProducer("VariableComputation",
                                RUNDATA = cms.InputTag("source", "RunData"),  
                                MWCHAMBERS = cms.InputTag("", ""),  
                                DWCTRACKS = cms.InputTag("", ""),                                  
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS"),
                                UserRecordCollectionName=cms.untracked.string("VariableUserRecords"),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                layerPositionFile=cms.string(layerPositionFile),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                NLayers=cms.untracked.int32(options.NLayers),
                                NColorsInputImage = cms.untracked.int32(-1),
                                CellEnergyCut = cms.untracked.double(0.5)
)


VariablesToPlot = ["xmean", "ymean", "NRechits", "E7_tot", "E19_tot", "EAll_tot", "EAllHG_tot", "EAllLG_tot", "EAllTOT_tot"]
VariablesToPlot += ["Ixx","Iyy","Izz","Ixy","Ixz","Iyz","depthX0","depthLambda0","showerStartDepth"]
for layer in range(1, 41):
    VariablesToPlot+=["EAll_layer%s"%layer, "EAll_layer%s"%layer, "EAll_layer%s"%layer]
process.observablentupler = cms.EDAnalyzer("NTupelizer",
                                USERRECORDS = cms.InputTag("variablecomputation","VariableUserRecords" ),
                                UserRecordKeys = cms.vstring(VariablesToPlot)
)


process.p = cms.Path( process.timingfilewriter * process.rawhitproducer * process.rechitproducer * process.rechitntupler * process.variablecomputation * process.observablentupler)
