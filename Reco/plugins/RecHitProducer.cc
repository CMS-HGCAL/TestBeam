#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
//#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBDigiCollections.h"

#define DEBUG


#ifdef DEBUG
#include <iostream>
#endif

class RecHitProducer : public edm::EDProducer{

      public:
             RecHitProducer(const edm::ParameterSet&);
             virtual void produce(edm::Event&, const edm::EventSetup&);             
      private:
              std::string outputCollectionName;     ///<label name of collection made by this producer
	edm::EDGetTokenT<HGCalTBDigiCollection> _digisToken;
};

RecHitProducer::RecHitProducer(const edm::ParameterSet& cfg)
	: outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
	  _digisToken(consumes<HGCalTBDigiCollection>(cfg.getParameter<edm::InputTag>("digiCollection")))
{
   produces <HGCalTBRecHitCollection>(outputCollectionName);

}

void RecHitProducer::produce(edm::Event& event, const edm::EventSetup&){
	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection); //auto_ptr are automatically deleted when out of scope
	
	edm::Handle<HGCalTBDigiCollection> digisHandle;
	event.getByToken(_digisToken, digisHandle);

	for(auto digi_itr = digisHandle->begin(); digi_itr != digisHandle->end(); ++digi_itr){
#ifdef DEBUG
		std::cout << *digi_itr << std::endl;
#endif
	}
    event.put(rechits, outputCollectionName);
  }
// Should there be a destructor ?? 
DEFINE_FWK_MODULE(RecHitProducer);
