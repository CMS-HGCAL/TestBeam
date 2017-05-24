#ifndef HGCALTBRECHITPRODUCER_H
#define HGCALTBRECHITPRODUCER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

//NUMBER_OF_SCA defined in HGCalTBSkiroc2CMS.h
struct pedestalChannel{
  HGCalTBDetId id;
  float pedHGMean[NUMBER_OF_SCA];
  float pedLGMean[NUMBER_OF_SCA];
  float pedHGRMS[NUMBER_OF_SCA];
  float pedLGRMS[NUMBER_OF_SCA];
};

class HGCalTBRawHitProducer : public edm::EDProducer
{

 public:
  HGCalTBRawHitProducer(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
 private:
  std::string m_outputCollectionName;
  std::string m_electronicMap;
  std::string m_pedestalLow_filename;
  std::string m_pedestalHigh_filename;

  std::map<int,pedestalChannel> m_pedMap; //key=10000*hexa+100*chip+chan
  
  edm::EDGetTokenT<HGCalTBSkiroc2CMSCollection> m_HGCalTBSkiroc2CMSCollection;

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  
};
#endif
