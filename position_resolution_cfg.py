import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys


options = VarParsing.VarParsing('standard')

####################################
# Options for reading in the data
options.register('fileDirectory',
                '/eos/user/t/tquast/data/Testbeam/July2017/DWC/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Path to the file from which the DWCs are read.'
                )

options.register('reportEvery',
                50000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.parseArguments()
input_data = []

#100 GeV pions
input_data.append((1159, "", 0, "1"))
input_data.append((1161, "", 0, "1"))
input_data.append((1162, "", 0, "1"))
input_data.append((1163, "", 0, "1"))
input_data.append((1165, "", 0, "1"))
input_data.append((1166, "", 0, "1"))   #check synchronisation of trigger
input_data.append((1171, "", 0, "1"))
input_data.append((1172, "", 0, "1"))
input_data.append((1173, "", 0, "1"))   #check synchronisation of trigger
input_data.append((1193, "", 0, "1"))



#150 GeV pions:
input_data.append((1051, "160717_1425", 1, "2"))
input_data.append((1067, "170717_0015", 0, "2"))
input_data.append((1069, "170717_0115", 0, "2"))
input_data.append((1070, "170717_0205", 0, "2"))
input_data.append((1071, "170717_0251", 0, "2"))
input_data.append((1072, "170717_0336", 0, "2"))
input_data.append((1073, "170717_0421", 0, "2"))
input_data.append((1074, "170717_0506", 0, "2"))
input_data.append((1076, "170717_0558", 0, "2"))
input_data.append((1079, "170717_0654", 0, "2"))
input_data.append((1082, "170717_0901", 1, "2"))
input_data.append((1088, "170717_1048", 1, "2"))
input_data.append((1090, "170717_1203", 0, "2"))
input_data.append((1096, "170717_1625", 0, "2"))

input_data.append((1134, "", 1, "2"))   #possibly remove last: careful check is necessary
input_data.append((1135, "", 0, "2"))   #check synchronisation of trigger
input_data.append((1136, "", 0, "2"))   #check synchronisation of trigger
input_data.append((1185, "", 0, "2"))   #check synchronisation of trigger
input_data.append((1192, "", 0, "2"))



#200 GeV pions
input_data.append((1132, "", 0, "3"))
input_data.append((1133, "", 0, "3"))   #check synchronisation of trigger
input_data.append((1190, "", 0, "3"))   #check synchronisation of trigger



#250 GeV pions
input_data.append((1189, "", 0, "4"))   #check synchronisation of trigger


#300 GeV pions
input_data.append((993, "150717_0800", 0, "5"))
input_data.append((994, "150717_0838", 0, "5"))
input_data.append((996, "150717_0924", 0, "5"))
input_data.append((997, "150717_1009", 0, "5"))
input_data.append((998, "150717_1100", 0, "5"))
input_data.append((999, "150717_1145", 0, "5"))
input_data.append((1000, "150717_1343", 0, "5"))
input_data.append((1001, "150717_1455", 0, "5"))
input_data.append((1002, "150717_1532", 0, "5"))
input_data.append((1003, "150717_1619", 0, "5"))
input_data.append((1004, "150717_1701", 0, "5"))
input_data.append((1012, "150717_1851", 0, "5"))
input_data.append((1013, "150717_1929", 0, "5"))
input_data.append((1025, "150717_2331", 0, "5"))
input_data.append((1027, "160717_0117", 0, "5"))
input_data.append((1028, "160717_0153", 0, "5"))
input_data.append((1029, "160717_0229", 0, "5"))
input_data.append((1030, "160717_0305", 0, "5"))
input_data.append((1032, "160717_0347", 0, "5"))
input_data.append((1033, "160717_0422", 0, "5"))
input_data.append((1034, "160717_0459", 0, "5"))
input_data.append((1035, "160717_0546", 0, "5"))
input_data.append((1036, "160717_0710", 0, "5"))
input_data.append((1039, "160717_0950", 0, "5"))
input_data.append((1040, "160717_1010", 1, "5"))
input_data.append((1045, "160717_1154", 0, "5"))
input_data.append((1050, "160717_1342", 1, "5"))


#150 GeV muons
#input_data.append((1140, "", 0, "12"))
#input_data.append((1181, "", 0, "12"))  #check synchronisation of trigger


#150 GeV electrons
input_data.append((1194, "", 0, "22"))


#250 GeV electrons
input_data.append((1038, "160717_0758", 0, "24"))



files = ["dwc_run_%s.root" % data[0] for data in input_data]
timingFileNames = ["RUN_%04d_%s_TIMING_RDOUT2.txt" % (data[0], data[1]) for data in input_data]


################################
# Setting an upper limit for the events to be processed, e.g. for debugging
options.maxEvents = -1
process = cms.Process("unpack")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

####################################
# Load the standard sequences
process.load('HGCal.StandardSequences.LocalReco_cff')
####################################


####################################
# Initialize the data read-in plugins
process.source = cms.Source("HGCalTBWireChamberSource",
    OutputCollectionName = cms.string("WireChambers"), 
    fileNames = cms.untracked.vstring(["file:%s/%s"%(options.fileDirectory, file) for file in files]),
    timingFileNames = cms.vstring(["/eos/user/t/tquast/data/Testbeam/July2017/Timing/%s"%file for file in timingFileNames]),
    skipFirstEventInDWCProducer = cms.vint32([data[2] for data in input_data]),
    runType = cms.vstring([data[3] for data in input_data]),
    performAlignment = cms.untracked.bool(True),
    alignmentParamaterFile = cms.untracked.string("/tmp/millepede.res") 

)


process.millepede_binarywriter.binaryFile = cms.string('/afs/cern.ch/user/t/tquast/millepede.bin')
process.millepede_binarywriter.nLayers = cms.int32(4)
process.millepede_binarywriter.MWCQualityCut = cms.bool(True)
process.millepede_binarywriter.makeTree = cms.untracked.bool(True)
process.millepede_binarywriter.MWCHAMBERS = cms.InputTag("source","WireChambers","unpack")
process.millepede_binarywriter.RUNDATA = cms.InputTag("source","RunData","unpack")
process.millepede_binarywriter.fittingMethod = cms.string("lineAnalytical")
process.millepede_binarywriter.binaryFile = cms.string("/tmp/millepede.bin")
                              


#tree file:
process.TFileService = cms.Service("TFileService", fileName = cms.string("outfile_DWCs.root"))

####################################
#add skip event exception which might occur for simulated samples because the last event is not properly passed forward
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.p = cms.Path(process.millepede_binarywriter*process.dwc_ntupelizer)

