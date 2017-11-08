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
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Reco/interface/PulseFitter.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iostream>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2F.h"
#include <sstream>


//#define DEBUG

class HGCalTBRecHitProducer : public edm::EDProducer
{
 public:
  HGCalTBRecHitProducer(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
 private:
  virtual void beginJob() override;
  std::string m_outputCollectionName;
  std::string m_electronicMap;
  bool m_maskNoisyChannels;
  std::string m_channelsToMask_filename;
  int m_NHexaBoards;
  double m_timeSample3ADCCut;
  
  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;
  std::vector<double> m_ADC_per_MIP;    
  std::vector<double> m_LG2HG_value;
  std::vector<double> m_TOT2LG_value;
  std::vector<double> m_highGainADCSaturation;
  std::vector<double> m_lowGainADCSaturation;

  bool performPulseFit;
  bool performAveraging;
  bool investigatePulseShape;

  std::map<int, TH2F*> shapesLG;
  std::map<int, TH2F*> shapesHG;
  std::map<int, TH2F*> ToARisevsTMaxLG;
  std::map<int, TH2F*> ToARisevsTMaxHG;
  std::map<int, TH2F*> ToAFallvsTMaxLG;
  std::map<int, TH2F*> ToAFallvsTMaxHG;
  std::map<int, TH2F*> TMaxHGvsTMaxLG;

  std::vector<int> m_noisyChannels;

  std::pair<double, double> CellCentreXY;
  HGCalTBCellVertices TheCell;

  #ifdef DEBUG
    int eventCounter;
  #endif

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  
};

#endif
