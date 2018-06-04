import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:


options.register('inputFiles',
                [''],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the input files.'
                )

options.register('processedFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/rechits/RECHITS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('outputFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where ntuples are stored')

options.register('physicsListUsed',
                "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Specify the used physics list to be passed forward to the run data object.'
                )

options.register('beamEnergy',
                250,
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

options.register('reportEvery',
                1000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

#hard coded configuration:
VariablesToPlot = ["eventID","run","pdgID","beamEnergy","configuration","runType","xmean","ymean","zmean","NRechits","NMIPHits","NNoisehits","25PercentQuantileRechitSpectrum","50PercentQuantileRechitSpectrum","75PercentQuantileRechitSpectrum","E1_tot","E7_tot","E19_tot","E37_tot","E61_tot","EAll_tot","EAllHG_tot","EAllLG_tot","EAllTOT_tot","depthX0"]

variables = []
for layer in range(1, 29):
    variables.append("NAll_layer%s"%layer)
    variables.append("EAll_layer%s"%layer)
    variables.append("E1PerE7_layer%s"%layer)
    variables.append("E7PerE19_layer%s"%layer)
    variables.append("E19PerE37_layer%s"%layer)
    variables.append("E37PerE61_layer%s"%layer)



electronicMap = "HGCal/CondObjects/data/map_CERN_Hexaboard_June2018_28Layers_dummy.txt"
setupConfiguration = 13
layerPositionFile = "/afs/cern.ch/user/t/tquast/CMSSW_9_3_0/src/HGCal/CondObjects/data/layer_distances_CERN_Hexaboard_June2018_28Layers_dummy.txt"  #attention: this path is hard coded and must be absolute
hgcalLayout = "HGCal/CondObjects/data/layerGeom_june2018_h2_28layers_dummy.txt"
areaSpecification = "H2"
NHexaBoards = 28
MaskNoisyChannels = False

options.parseArguments()
print options


################################
process = cms.Process("gensim")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

################################
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################


process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.processedFile),                                  
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*',
                                                                         'keep *_*_HGCalTBDATURATracks_*',
                                                                         'keep *_*_FullRunData_*')
)

process.source = cms.Source("HGCalTBGenSimSource",
                        fileNames=cms.untracked.vstring(["file:%s" % file for file in options.inputFiles]),
                        RechitOutputCollectionName = cms.string('HGCALTBRECHITS'), 
                        produceDATURATracksInsteadOfDWCs = cms.untracked.bool(True),
                        DWCOutputCollectionName = cms.string(''), 
                        DATURAOutputCollectionName = cms.string('HGCalTBDATURATracks'), 
                        RunDataOutputCollectionName = cms.string('FullRunData'), 
                        e_mapFile_CERN = cms.untracked.string(electronicMap),
                        layerPositionFile=cms.string(layerPositionFile),
                        MaskNoisyChannels=cms.untracked.bool(MaskNoisyChannels),
                        ChannelsToMaskFileName=cms.untracked.string(""),
                        beamEnergy=cms.untracked.double(options.beamEnergy),
                        beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),                        
                        energyNoise=cms.untracked.double(0.0),  #indicated in MIPs
                        setupConfiguration=cms.untracked.uint32(setupConfiguration),
                        energyNoiseResolution=cms.untracked.double(1./6.), #indicated in MIPs
                        GeVToMip=cms.untracked.double(1./(84.9*pow(10.,-6))),       #assume one MIP = 84.9keV
                        areaSpecification = cms.untracked.string(areaSpecification),
                        physicsListUsed = cms.untracked.string(options.physicsListUsed),
                        wc_resolutions = cms.untracked.vdouble(4*[0.2])        #set to the expected resolutions according to the manual
                        )


rundata_tag = cms.InputTag("source", "FullRunData" )
rechit_tag = cms.InputTag("source","HGCALTBRECHITS" )

process.rechitntupler = cms.EDAnalyzer("RecHitNtupler",
                                       InputCollection=rechit_tag,
                                       RUNDATA = rundata_tag,
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       layerPositionFile = cms.untracked.string(layerPositionFile),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(True),
                                       MipThreshold=cms.untracked.double(5.0),
                                       NoiseThreshold=cms.untracked.double(0.5)
)

process.variablecomputation = cms.EDProducer("VariableComputation",
                                RUNDATA = rundata_tag,  
                                MWCHAMBERS = cms.InputTag("","" ), 
                                DWCTRACKS = cms.InputTag("","" ),                                 
                                HGCALTBRECHITS = rechit_tag,
                                UserRecordCollectionName=cms.untracked.string("VariableUserRecords"),
                                ElectronicMap = cms.untracked.string(electronicMap),
                                DetectorLayout=cms.untracked.string(hgcalLayout),
                                layerPositionFile=cms.string(layerPositionFile),
                                NHexaBoards=cms.untracked.int32(NHexaBoards),
                                NLayers=cms.untracked.int32(NHexaBoards),
                                NColorsInputImage = cms.untracked.int32(-1),
                                CellEnergyCut = cms.untracked.double(0.5)
)

process.ntupelizer = cms.EDAnalyzer("NTupelizer",
                                USERRECORDS = cms.InputTag("variablecomputation","VariableUserRecords" ),
                                UserRecordKeys = cms.vstring(VariablesToPlot)
)


####################################

process.p = cms.Path( process.rechitntupler * process.variablecomputation * process.ntupelizer)


process.end = cms.EndPath(process.output)
