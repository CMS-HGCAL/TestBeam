#ifndef HGCALTBRECHITPRODUCER_H
#define HGCALTBRECHITPRODUCER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/CondObjects/interface/HGCalTBADCConversionsMap.h"
#include "HGCal/CondObjects/interface/HGCalTBADCConversionsMap_perChannel.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2F.h"


enum PreselectionMethod {
  TB2017 = 0,
  TB2018, 
  NONE
};

class HGCalTBRecHitProducer : public edm::EDProducer
{
 public:
  HGCalTBRecHitProducer(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
 private:
  virtual void beginJob() override;
  std::string m_outputCollectionName;
  std::string m_electronicMap;
  std::string m_detectorLayoutFile;
  std::string m_adcCalibrationsFile;
  bool m_calibrationPerChannel;
  int m_expectedMaxTimeSample;
  double m_maxADCCut;
  bool m_useCalibration;
  bool m_subtractCommonMode;
  int m_TSForCommonModeNoiseSubtraction;

  std::map<int, TH2F*> shapesLG;
  std::map<int, TH2F*> shapesHG;
  
  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  float m_fittingTime;

  std::pair<Float16_t, Float16_t> CellCentreXY;
  HGCalTBCellVertices TheCell;

  std::string m_preselectionMethod;
  PreselectionMethod _preselectionMethod;

  struct {
    HGCalElectronicsMap emap_;
    HGCalTBDetectorLayout layout_;
    HGCalTBADCConversionsMap adccalibmap_;
    HGCalTBADCConversionsMap_perChannel adccalibmap_perchannel_;
  } essource_;
  
};

#endif