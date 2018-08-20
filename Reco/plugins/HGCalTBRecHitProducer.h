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
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectROOTIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"

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
  int m_expectedMaxTimeSample;
  double m_maxADCCut;
  bool m_useCalibration;
  
  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  float m_fittingTime;

  std::pair<double, double> CellCentreXY;
  HGCalTBCellVertices TheCell;

  HGCalElectronicsMap m_emap;
  HGCalTBDetectorLayout m_layout;
  HGCalTBADCConversionsMap m_adccalibmap;

};

#endif
