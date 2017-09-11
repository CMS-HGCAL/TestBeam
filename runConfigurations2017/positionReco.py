import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'


options.register('dataFile',
                 '/home/tquast/tb2017/reconstructedFiles/reco_1305.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'folder containing raw input')

options.register('outputFile',
                 '/home/tquast/tb2017/analysis/positionResolution_1305.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output folder where analysis output are stored')

options.register('reportEvery',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )


####################################
# Options related to the experimental setup
options.register('configuration',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 '1 is the July setup'
                 )


options.register('alignmentParameterFiles',
                 '',
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Alignment parameter files as obtained from the (mille)pede framework.'
                )


#Specific options for the extra(inter)polation of impact positions in the MillepedeBinaryWrite and positionResolutionAnalyzer
options.register('useMWCReference',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 'Are the multi wire chamber information used for the exrapolation of referemce impact points?')

options.register('fittingMethod',
                'lineAnalytical',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Model of the electron tracks.'    
    )


####################################
#Specific options for reconstruction of impact positions in the MillepedeBinaryWrite and positionResolutionAnalyzer
options.register('considerationMethod',
                 'closest19',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Possible arguments are: all, closest7, closest19, clusters.'
                )

options.register('weightingMethod',
                 'linearWeighting',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Possible arguments are: squaredWeighting, linearWeighting, logWeighting_<a>_1.0, .... '
                )




options.maxEvents = -1


options.parseArguments()
print options



################################
process = cms.Process("positionresolutionanalysis")
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

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))

####################################
process.position_resolution_analyzer = cms.EDAnalyzer("Position_Resolution_Analyzer",
                                alignmentParameterFiles = cms.vstring(options.alignmentParameterFiles),
                                considerationMethod = cms.string(options.considerationMethod),
                                weightingMethod = cms.string(options.weightingMethod),
                                fittingMethod = cms.string(options.fittingMethod),
                                layers_config  = cms.int32(1),
                                ADC_per_MIP = cms.vdouble([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]),
                                nLayers = cms.int32(6),
                                SensorSize = cms.int32(133),
                                useMWCReference = cms.bool(options.useMWCReference),
                                RUNDATA = cms.InputTag("source","RunData" ), 
                                MWCHAMBERS = cms.InputTag("wirechamberproducer","DelayWireChambers" ), 
                                HGCALTBRECHITS = cms.InputTag("rechitproducer","HGCALTBRECHITS" ),
                              )


####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
process.load('HGCal.StandardSequences.RawToDigi_cff')
####################################


process.p = cms.Path(process.position_resolution_analyzer)

