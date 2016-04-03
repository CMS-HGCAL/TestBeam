import FWCore.ParameterSet.Config as cms

pedestals = cms.EDAnalyzer('Pedestals',
                           digiCollection = cms.InputTag('hgcaltbdigis'),
)
