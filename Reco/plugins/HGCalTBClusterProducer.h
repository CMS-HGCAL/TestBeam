#ifndef HGCALTBCLUSTERPRODUCER_H
#define HGCALTBCLUSTERPRODUCER_H
/** \class Reco/plugins/HGCalTBClusterProducer.h 
    \brief

    \author Arnaud Steen
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"

#include <iostream>

class HGCalTBClusterProducer : public edm::EDProducer
{

 public:
  HGCalTBClusterProducer(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
 private:
  virtual void beginJob() override;
  std::string m_electronicMap;
  std::string m_detectorLayoutFile;
  int m_sensorSize;
  std::string m_outputCollectionName;
  std::string m_outputCollectionName7;
  std::string m_outputCollectionName19;
  bool m_runDynamicCluster;
  bool m_runCluster7;
  bool m_runCluster19;
  double m_minEnergy;
	
  edm::EDGetTokenT<HGCalTBRecHitCollection> m_HGCalTBRecHitCollection;

  struct {
    HGCalElectronicsMap emap;
    HGCalTBDetectorLayout layout;
  } m_essource;

  void buildCluster(HGCalTBRecHitCollection rechits, std::vector<HGCalTBDetId> &temp, std::vector<HGCalTBDetId> &clusterDetIDs, int maxDistance);
  void createDynamicClusters(HGCalTBRecHitCollection rechits, std::vector<reco::HGCalTBCluster> &clusterCol);
  void createSeededClusters(HGCalTBRecHitCollection rechits, reco::HGCalTBCluster &cluster, int maxDistance);
};

#endif
