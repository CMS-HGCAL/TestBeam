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

#include <iostream>

class HGCalTBClusterProducer : public edm::EDProducer
{

public:
	HGCalTBClusterProducer(const edm::ParameterSet&);
	virtual void produce(edm::Event&, const edm::EventSetup&);
private:
	std::string _elecMapFile;
	std::string _outputCollectionName;
	std::string _outputCollectionName7;
	std::string _outputCollectionName19;
	edm::EDGetTokenT<HGCalTBRecHitCollection> _rechitToken;
	HGCalElectronicsMap _elecMap;
	
	bool _runDynamicCluster;
	bool _runCluster7;
	bool _runCluster19;
	int _sensorSize;
	int _maxTransverse;
	double _minEnergy;

	std::vector<double> _layerZPositions;

	void buildCluster(HGCalTBRecHitCollection rechits, std::vector<HGCalTBDetId> &temp, std::vector<HGCalTBDetId> &clusterDetIDs);
	void createDynamicClusters(HGCalTBRecHitCollection rechits, std::vector<reco::HGCalTBCluster> &clusterCol);
	void createSeededClusters(HGCalTBRecHitCollection rechits, reco::HGCalTBCluster &cluster, int maxDist);
};

#endif
