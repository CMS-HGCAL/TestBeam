import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbrechitproducer_cfi import *

LocalRecoSeq  = cms.Sequence(hgcaltbrechits)

