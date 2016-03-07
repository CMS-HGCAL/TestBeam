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

class HGCalTBRecHitProducer : public edm::EDProducer{

      public:
             HGCalTBRecHitProducer(const edm::ParameterSet&);
             virtual void produce(edm::Event&, const edm::EventSetup&);             
      private:
              std::string outputCollectionName;     ///<label name of collection made by this producer
	edm::EDGetTokenT<HGCalTBDigiCollection> _digisToken;
};

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg)
	: outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
	  _digisToken(consumes<HGCalTBDigiCollection>(cfg.getParameter<edm::InputTag>("digiCollection")))
{
   produces <HGCalTBRecHitCollection>(outputCollectionName);

}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup&){

	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection); //auto_ptr are automatically deleted when out of scope
	
	edm::Handle<HGCalTBDigiCollection> digisHandle;
	event.getByToken(_digisToken, digisHandle);

	for(auto digi_itr = digisHandle->begin(); digi_itr != digisHandle->end(); ++digi_itr){
#ifdef DEBUG
		std::cout << "[RECHIT]" << *digi_itr << std::endl;
#endif
		SKIROC2DataFrame digi(*digi_itr);
		unsigned int nSamples = digi.samples();
		for(unsigned int iSample = 0; iSample < nSamples; ++iSample){
			HGCalTBRecHit recHit(digi.detid(), digi[iSample].adc(), digi[iSample].tdc()); ///\todo use calibration!
			rechits->push_back(recHit);
		}
	}
    event.put(rechits, outputCollectionName);
  }
// Should there be a destructor ?? 
DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
