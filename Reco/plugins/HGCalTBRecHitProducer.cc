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

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup){

	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection); //auto_ptr are automatically deleted when out of scope
	
	edm::Handle<HGCalTBDigiCollection> digisHandle;
	event.getByToken(_digisToken, digisHandle);

#ifdef ESPRODUCER
// the conditions should be defined by ESProduers and available in the iSetup
	edm::ESHandle<HGCalCondPedestals> pedestalsHandle;
	//iSetup.get<HGCalCondPedestals>().get(pedestalsHandle); ///\todo should we support such a syntax?
	HGCalCondPedestals pedestals = iSetup.get<HGCalCondPedestals>();

	edm::ESHandle<HGCalCondGains>     adcToGeVHandle;
	//iSetup.get<HGCalCondGains>().get(adcToGeVHandle);
	HGCalCondGains adcToGeV = iSetup.get<HGCalCondGains>();
#else
	HGCalCondGains adcToGeV;
	HGCalCondPedestals pedestals;
	HGCalCondObjectTextIO condIO(HGCalTBNumberingScheme::scheme());
    //HGCalElectronicsMap emap;
//    assert(io.load("mapfile.txt",emap)); ///\todo to be trasformed into exception
	assert(condIO.load("pedestals.txt", pedestals));
	assert(condIO.load("gains.txt", adcToGeV));	
	///\todo check if reading the conditions from file some channels are not in the file!	  
#endif


	for(auto digi_itr = digisHandle->begin(); digi_itr != digisHandle->end(); ++digi_itr){
#ifdef DEBUG
		std::cout << "[RECHIT PRODUCER: digi]" << *digi_itr << std::endl;
#endif

//------------------------------ this part should go into HGCalTBRecoAlgo class

		SKIROC2DataFrame digi(*digi_itr);
		unsigned int nSamples = digi.samples();
		for(unsigned int iSample = 0; iSample < nSamples; ++iSample){
			
			float energy = (digi[iSample].adc() - pedestals.get(digi.detid())->value) * adcToGeV.get(digi.detid())->value;

			HGCalTBRecHit recHit(digi.detid(), energy, digi[iSample].tdc()); ///\todo use time calibration!

#ifdef DEBUG
			std::cout << recHit << std::endl;
#endif
			rechits->push_back(recHit);
		}
	}
    event.put(rechits, outputCollectionName);
  }
// Should there be a destructor ?? 
DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
