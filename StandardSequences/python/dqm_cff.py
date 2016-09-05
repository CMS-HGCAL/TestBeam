import FWCore.ParameterSet.Config as cms

from HGCal.RawToDigi.hgcaltbdigisplotter_cfi import *
from HGCal.Reco.hgcaltbrechitplotter_cfi import *

DQMDigiSeq   = cms.Sequence(hgcaltbdigisplotter)
DQMRecHitSeq =  cms.Sequence(hgcaltbrechitsplotter)

