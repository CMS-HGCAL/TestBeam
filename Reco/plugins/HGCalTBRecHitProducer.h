#ifndef HGCALTBRECHITPRODUCER_H
#define HGCALTBRECHITPRODUCER_H
/** \class Reco/plugins/HGCalTBRecHitProducer.h HGCalTBRecHitProducer HGCalTBRecHitProducer
	\brief

	\author Shervin Nourbakhsh
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
//#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBDigiCollections.h"


// here are defined objects containing pedestals and ADCtoGeV factors
#include "HGCal/CondObjects/interface/HGCalCondObjects.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalTBNumberingScheme.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
//#define DEBUG


#ifdef DEBUG
#include <iostream>
#endif

class HGCalTBRecHitProducer : public edm::EDProducer
{

public:
	HGCalTBRecHitProducer(const edm::ParameterSet&);
	virtual void produce(edm::Event&, const edm::EventSetup&);
private:
	std::string outputCollectionName;     ///<label name of collection made by this producer
	edm::EDGetTokenT<HGCalTBDigiCollection> _digisToken;
	std::string _pedestalLow_filename, _pedestalHigh_filename;
	std::string _gainsLow_filename, _gainsHigh_filename;
	int _adcSaturation;
	std::vector<double> _LG2HG_value;
	std::string _mapFile;
	int _layers_config;
	struct {
	  HGCalElectronicsMap emap_;
        } essource_;
};



#endif
