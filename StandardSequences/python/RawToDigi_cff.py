import FWCore.ParameterSet.Config as cms

from HGCal.RawToDigi.hgcaltbdigis_cfi import *
from HGCal.RawToDigi.hgcaltbtracks_cfi import *

RawToDigiSeq = cms.Sequence(hgcaltbdigis * hgcaltbtracks)
