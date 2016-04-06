#include "HGCal/Reco/plugins/HGCalTBRecHitProducerLowGain.h"
#include <iostream>
using namespace std;

HGCalTBRecHitProducerLowGain::HGCalTBRecHitProducerLowGain(const edm::ParameterSet& cfg)
	: outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
	  _digisToken(consumes<HGCalTBDigiCollection>(cfg.getParameter<edm::InputTag>("digiCollection")))
{
	produces <HGCalTBRecHitCollection>(outputCollectionName);

}

void HGCalTBRecHitProducerLowGain::produce(edm::Event& event, const edm::EventSetup& iSetup)
{

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
	assert(condIO.load("Ped_LowGain_Test_8272.txt", pedestals));
	assert(condIO.load("Gain_Test.txt", adcToGeV));
	///\todo check if reading the conditions from file some channels are not in the file!
#endif


	for(auto digi_itr = digisHandle->begin(); digi_itr != digisHandle->end(); ++digi_itr) {
#ifdef DEBUG
		std::cout << "[RECHIT PRODUCER: digi]" << *digi_itr << std::endl;
#endif

//------------------------------ this part should go into HGCalTBRecoAlgo class

		SKIROC2DataFrame digi(*digi_itr);
		unsigned int nSamples = 1;

		// if there are more than 1 sample, we need to define a reconstruction algorithm
		// now taking the first sample
		for(unsigned int iSample = 0; iSample < nSamples; ++iSample) {

//			float energy = (digi[iSample].adc() - pedestals.get(digi.detid())->value) * adcToGeV.get(digi.detid())->value;
//                        float energy = (digi[iSample].adc() - pedestals.get(digi.detid())->value);
                        float energy = (digi[iSample].adc() - pedestals.get(digi.detid())->value);
			HGCalTBRecHit recHit(digi.detid(), energy, 0.); ///\todo use time calibration!

#ifdef DEBUG
			std::cout << recHit << std::endl;
#endif
			if(iSample == 0) rechits->push_back(recHit); ///\todo define an algorithm for the energy if more than 1 sample, code inefficient
		}
	}
	event.put(rechits, outputCollectionName);
}
// Should there be a destructor ??
DEFINE_FWK_MODULE(HGCalTBRecHitProducerLowGain);
