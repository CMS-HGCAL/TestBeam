import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbshower_cfi import *

ShowerRecoSeq  = cms.Sequence(hgcaltbshower)

