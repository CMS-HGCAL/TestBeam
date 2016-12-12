import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbtracking_cfi import *

TrackingRecoSeq  = cms.Sequence(hgcaltbtrackingexample)
TrackAnaRecoSeq  = cms.Sequence(hgcaltbtrackanalyzer)
