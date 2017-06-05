#ifndef HGCALTBRECHITPRODUCER_H
#define HGCALTBRECHITPRODUCER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

class HGCalTBRecHitProducer : public edm::EDProducer
{
 public:
  HGCalTBRecHitProducer(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
 private:
  virtual void beginJob() override;
  std::string m_pedestalHigh_filename;
  std::string m_pedestalLow_filename;
  std::string m_outputCollectionName;
  std::string m_electronicMap;
  double m_commonModeThreshold;
  int m_highGainADCSaturation;
  int m_lowGainADCSaturation;
  bool m_keepOnlyTimeSample3;

  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;
  std::vector<double> m_LG2HG_value;
  std::vector<double> m_TOT2LG_value;

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  struct pedestalChannel{
    HGCalTBDetId id;
    float pedHGMean[NUMBER_OF_TIME_SAMPLES];
    float pedLGMean[NUMBER_OF_TIME_SAMPLES];
    float pedHGRMS[NUMBER_OF_TIME_SAMPLES];
    float pedLGRMS[NUMBER_OF_TIME_SAMPLES];
  };
  std::map<int,pedestalChannel> m_pedMap; //key=10000*hexa+100*chip+chan
  
  struct commonModeNoise{
  commonModeNoise():fullHG(0),halfHG(0),mouseBiteHG(0),outerHG(0),fullLG(0),halfLG(0),mouseBiteLG(0),outerLG(0),fullCounter(0),halfCounter(0),mouseBiteCounter(0),outerCounter(0){;}
    float fullHG,halfHG,mouseBiteHG,outerHG;
    float fullLG,halfLG,mouseBiteLG,outerLG;
    int fullCounter,halfCounter,mouseBiteCounter,outerCounter;
  };

};

#endif