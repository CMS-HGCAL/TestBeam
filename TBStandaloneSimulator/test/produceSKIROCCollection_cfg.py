# ---------------------------------------------------------------------------
# Read sim hits from sim root file and store in them in edm::Events 
# together with SKIROC2DataFrames
# Created April 2016 HBP 
# ---------------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms

processName = "HGC"

process = cms.Process(processName)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# -------------------------------------------------
# read file names from filelist and noise_filelist
# -------------------------------------------------
from string import strip
import os, sys

filelist = map(lambda x: "file:%s" % x, 
               filter(lambda x: len(x) > 0 and x[0] != "#",
                      map(strip, open("filelist").readlines())))

# read file names from the noise filelist
if not os.path.exists("noise_filelist"):
    print "** noise_filelist NOT FOUND - creating an empty file"
    os.system('echo "#" > noise_filelist')
noisefilelist = map(lambda x: "file:%s" % x, 
                    filter(lambda x: len(x) > 0 and x[0] != "#",
                           map(strip, open("noise_filelist").readlines())))
# -------------------------------------------------
process.source = cms.Source ("HGCSimDigiSource",
                             runNumber  = cms.untracked.int32(101),
                             maxEvents  = cms.untracked.int32(-1),
                             minADCCount= cms.untracked.int32(1),
                             # this value will actually give the energy in keV
                             ADCperMeV  = cms.untracked.double(1000),
                             fileNames  = cms.untracked.vstring(filelist),
                             noiseFileNames = 
                             cms.untracked.vstring(noisefilelist)
                             )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string
                               ("HGC_electrons_8GeV_2016_4layer_20k_sim.root")
                               )

process.outpath = cms.EndPath(process.out)
