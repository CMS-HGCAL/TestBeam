import FWCore.ParameterSet.Config as cms

imaging_tuple_writer = cms.EDAnalyzer("Imaging_Tuple_Writer",
                                RUNDATA = cms.InputTag("source","RunData","unpack" ), 
                                HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                ADC_per_MIP = cms.vdouble([17.31, 17.12, 16.37, 17.45, 17.31, 16.98, 16.45, 16.19, 17.55, 17.19, 16.99, 17.92, 15.95, 16.64, 16.79, 15.66]),
                                e_mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt')
                              )
