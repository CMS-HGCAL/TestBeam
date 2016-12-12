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

#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrackCollection.h"

#include <iostream>

class HGCalTBCaloTrackProducer : public edm::EDProducer
{

public:
	HGCalTBCaloTrackProducer(const edm::ParameterSet&);
	virtual void produce(edm::Event&, const edm::EventSetup&);
private:
	std::string _outputCollectionName;
	edm::EDGetToken HGCalTBClusterCollection_;
	
	bool _doTrackCleaning;
	double _maxDistanceToRecoTrack;
	int _minTouchedLayers;
	double _minEnergy;
	double _maxEnergy;

	std::vector<double> _layerZPositions;
};

#endif
