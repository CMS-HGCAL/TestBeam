import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbtracking_cfi import *

TrackAnaRecoSeq  = cms.Sequence(hgcaltbtrackanalyzer)
