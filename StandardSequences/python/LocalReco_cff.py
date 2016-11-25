import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbrechitproducer_cfi import *
from HGCal.Reco.hgcaltbclusterproducer_cfi import *

LocalRecoSeq  = cms.Sequence(hgcaltbrechits)
ClusterRecoSeq  = cms.Sequence(hgcaltbclusters)

