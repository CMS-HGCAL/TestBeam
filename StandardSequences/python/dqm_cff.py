import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbrechitplotter_cfi import *

DQMRecHitSeq =  cms.Sequence(hgcaltbrechitsplotter)

