import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('runNumber',
                496,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'RunNumber.'
                )

options.register('beamEnergy',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Beam energy.'
                )

options.register('beamParticlePDGID',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Beam particles PDG ID.'
                )

options.register('runType',
                 "Beam",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Run type: Pedestal, Beam, Simulation.'
                )

options.register('setupConfiguration',
                5,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'setupConfiguration (1: July - 4: 20 Layers in October in H6A".'
                )

options.register('NHexaBoards',
                4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('SubtractPedestal',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Subtract the pedestals.'
                )

options.register('electronicMap',
                 'map_CERN_Hexaboard_July_6Layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('adcCalibrations',
                 'hgcal_calibration.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal ADC to MIP calibration file in HGCal/CondObjects/data/')

options.register('hgcalLayout',
                 'layerGeom_oct2017_h2_17layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')


options.register('layerPositionFile',
                 '/afs/cern.ch/user/t/tquast/CMSSW_8_3_0/src/HGCal/CondObjects/data/layer_distances_DESY_March2018_config4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('pathsToMIPWindowFiles',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the file containing the corresponding cuts on the reconstructed impact positions on DWC E.'
                )

options.register('NBins',
                60,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Bins for the DWC search window.'
                )

options.register('DimXDUT',
                30.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'MIP search window on DWCE in x-coordinate [mm].'
                )

options.register('DimYDUT',
                30.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'MIP search window on DWCE in y-coordinate [mm].'
                )

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
process = cms.Process("mipfindinganalysis")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################

#only read time stamps for beam runs
readTimeStamps = False
if options.beamParticlePDGID==0:
    readTimeStamps = False




process.source = cms.Source("HGCalTBRawDataSource",
                            ElectronicMap=cms.untracked.string(electronicMap),
                            fileNames=cms.untracked.vstring("file:/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/ORM_raw/HexaData_Run%04d.raw"%options.runNumber),
                            OutputCollectionName=cms.untracked.string("skiroc2cmsdata"),
                            NumberOf32BitsWordsPerReadOut=cms.untracked.uint32(30787),
                            NumberOfBytesForTheHeader=cms.untracked.uint32(12),          #for the new headers/trailers from run 1241 onward: 12
                            NumberOfBytesForTheTrailer=cms.untracked.uint32(4),         #for the new headers/trailers from run 1241 onward: 4
                            NumberOfBytesForTheEventTrailers=cms.untracked.uint32(12),   #for the new headers/trailers from run 1241 onward: 12
                            NSkipEvents=cms.untracked.uint32(0),
                            ReadTimeStamps=cms.untracked.bool(readTimeStamps),
                            DataFormats=cms.untracked.uint32(1),
                            timingFiles=cms.vstring(),
                            beamEnergy=cms.untracked.double(options.beamEnergy),
                            beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),
                            runType=cms.untracked.string(options.runType),
                            setupConfiguration=cms.untracked.uint32(options.setupConfiguration)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/quickAnalysis_thorben/analysis/analysed_%04d.root"%options.runNumber))


process.rawhitproducer = cms.EDProducer("HGCalTBRawHitProducer",
                                        InputCollection=cms.InputTag("source","skiroc2cmsdata"),
                                        OutputCollectionName=cms.string("HGCALTBRAWHITS"),
                                        GlobalTimestampCollectionName=cms.string("HGCALGLOBALTIMESTAMPS"),
                                        ElectronicMap=cms.untracked.string(electronicMap),
                                        SubtractPedestal=cms.untracked.bool(bool(1)),
                                        MaskNoisyChannels=cms.untracked.bool(bool(0)),
                                        HighGainPedestalFileName=cms.untracked.string("/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/quickAnalysis_thorben/pedestals/pedestalHG_%04d.txt"%(options.runNumber)),
                                        LowGainPedestalFileName=cms.untracked.string("/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/quickAnalysis_thorben/pedestals/pedestalLG_%04d.txt"%(options.runNumber)),
                                        NoisyChannelsFileName=cms.untracked.string("/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/quickAnalysis_thorben/pedestals/noisyChannels_%04d.txt"%(options.runNumber))
)

process.rawhitplotter = cms.EDAnalyzer("RawHitPlotter",
                                       InputCollection=cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(False),
                                       SubtractCommonMode=cms.untracked.bool(True)
)

process.rechitproducer = cms.EDProducer("HGCalTBRecHitProducer",
                                        OutputCollectionName = cms.string('HGCALTBRECHITS'),
                                        InputCollection = cms.InputTag("rawhitproducer","HGCALTBRAWHITS"),
                                        RUNDATA = cms.InputTag("source", "RunData" ), 
                                        GlobalTimestampCollectionName=cms.InputTag("rawhitproducer","HGCALGLOBALTIMESTAMPS"),
                                        ElectronicsMap = cms.untracked.string(electronicMap),
                                        DetectorLayout = cms.untracked.string(hgcalLayout),
                                        ADCCalibrations = cms.untracked.string(adcCalibrations),                                       
                                        MaskNoisyChannels=cms.untracked.bool(bool(0)),
                                        ChannelsToMaskFileName=cms.untracked.string("/eos/cms/store/group/dpg_hgcal/tb_hgcal/desy_march2018/quickAnalysis_thorben/pedestals/noisyChannels_%04d.txt"%(options.runNumber)),
                                        NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                        TimeSample3ADCCut = cms.untracked.double(15.),
                                        investigatePulseShape = cms.untracked.bool(True),
                                        timingNetworks = cms.untracked.string("")
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


process.mipfindinganalysis = cms.EDAnalyzer("MIPFinder",
                                RUNDATA = cms.InputTag("source", "RunData" ), 
                                MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                DWCTRACKS = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                HGCALTBCOMMONMODENOISE = cms.InputTag("rechitproducer","HGCALTBCOMMONMODENOISEMAP" ),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                n_bins_DWCE = cms.int32(options.NBins),
                                max_dim_x_DUT = cms.double(options.DimXDUT),
                                max_dim_y_DUT = cms.double(options.DimYDUT),
                                pathsToMIPWindowFiles = cms.vstring(options.pathsToMIPWindowFiles),
                                commonModeNoiseRejectionType = cms.int32(0)       #0: none, else 1-..., default: 0 

                              )

process.variablecomputation = cms.EDProducer("VariableComputation",
                                RUNDATA = cms.InputTag("source", "RunData" ), 
                                MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                DWCTRACKS = cms.InputTag("dwctrackproducer","HGCalTBDWCTracks" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                                UserRecordCollectionName=cms.untracked.string("VariableUserRecords"),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                layerPositionFile=cms.string(options.layerPositionFile),
                                NHexaBoards=cms.untracked.int32(options.NHexaBoards),
                                NLayers=cms.untracked.int32(options.NHexaBoards),
                                DNNInputFile = cms.untracked.string(""),
                                NColorsInputImage = cms.untracked.int32(-1)
                              )

process.ntupelizer = cms.EDAnalyzer("NTupelizer",
                                USERRECORDS = cms.InputTag("variablecomputation","VariableUserRecords" ),
                                UserRecordKeys = cms.vstring(["NRechits", "I_EV1", "I_EV2", "E1_tot", "E7_tot", "E19_tot", "E19_tot", "E37_tot", "EAll_tot", "depthX0", "NAll_layer1", "NAll_layer2", "NAll_layer3", "E1PerE7_layer1", "E1PerE7_layer2", "E1PerE7_layer3", "E7PerE19_layer1", "E7PerE19_layer2", "E7PerE19_layer3"])
                            )



####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.rawhitproducer*process.rawhitplotter*process.rechitproducer*process.rechitplotter*process.mipfindinganalysis*process.variablecomputation*process.ntupelizer)

