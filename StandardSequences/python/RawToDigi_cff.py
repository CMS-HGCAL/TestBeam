import FWCore.ParameterSet.Config as cms

from HGCal.RawToDigi.hgcaltbdigis_cfi import *

RawToDigiSeq = cms.Sequence(hgcaltbdigis)
