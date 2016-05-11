// -*- C++ -*-
//
// Package:    HGCal/Calibration
// Class:      Pedestals
//
/**\class Pedestals Pedestals.cc HGCal/Pedestals/plugins/Pedestals.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shervin Nourbakhsh
//         Created:  Sun, 03 Apr 2016 09:12:17 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "HGCal/DataFormats/interface/HGCalTBDigiCollections.h"
#include "HGCal/CondObjects/interface/HGCalCondObjects.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalTBNumberingScheme.h"

#include <cassert>

//
// class declaration
//
/// \todo calculated and define pedestals for all gains

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Pedestals : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
	explicit Pedestals(const edm::ParameterSet&);
	~Pedestals();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

	// ----------member data ---------------------------

	edm::EDGetTokenT<HGCalTBDigiCollection> _digisToken;

	typedef struct {
		float sum;
		float sum2;
		unsigned int n;
	} pedestalSum_t;

	typedef std::map< HGCalTBDetId, pedestalSum_t > pedestalMap_t;
	pedestalMap_t pedestals;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Pedestals::Pedestals(const edm::ParameterSet& iConfig) :
	_digisToken(consumes<HGCalTBDigiCollection>(iConfig.getParameter<edm::InputTag>("digiCollection")))
{
	//now do what ever initialization is needed
	usesResource("TFileService");

}


Pedestals::~Pedestals()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Pedestals::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	edm::Handle<HGCalTBDigiCollection> digisHandle;
	iEvent.getByToken(_digisToken, digisHandle);

	HGCalCondPedestals pedestals_cond(HGCalTBNumberingScheme::scheme(), 0);
	HGCalCondNoise noise_cond(HGCalTBNumberingScheme::scheme(), 0);

	HGCalCondObjectTextIO condIO(HGCalTBNumberingScheme::scheme());
	HGCalElectronicsMap emap;

	condIO.load("CondObjects/data/map_FNAL_2.txt", emap);

	for(auto digi_itr = digisHandle->begin(); digi_itr != digisHandle->end(); ++digi_itr) {

		SKIROC2DataFrame digi(*digi_itr);
		HGCalTBDetId detId = digi.detid();

		unsigned int nSamples = digi.samples();

		/**** WARNING \todo multiple samples are not yet allowed */
		nSamples = 1;
		/**** */

		pedestalMap_t::iterator thisPedestal_itr = pedestals.find(detId);
		if(thisPedestal_itr == pedestals.end()) {
			pedestalSum_t empty_ped;
			empty_ped.sum = 0.;
			empty_ped.n = 0;

			pedestals[detId] = empty_ped;
			thisPedestal_itr = pedestals.find(detId);
			assert(thisPedestal_itr != pedestals.end());
		}

		// now sum the values from different events
		for(unsigned int iSample = 0; iSample < nSamples; ++iSample) {
			thisPedestal_itr->second.sum += digi[iSample].adcLow();
			thisPedestal_itr->second.sum2 += digi[iSample].adcLow() * digi[iSample].adcLow();
			++(thisPedestal_itr->second.n);
#ifdef DEBUG
			std::cout << "[PEDESTAL PRODUCER: digi]" << *digi_itr << std::endl;
			std::cout << "                          ** " << thisPedestal_itr->second.sum << std::endl;
#endif

		}
	}

	// now I have the pedestal average, I have to print them
	for(const auto pedestal : pedestals) {
		float pedestal_mean = pedestal.second.sum / pedestal.second.n;
		pedestals_cond.set(pedestal.first, pedestal_mean);
		noise_cond.set(pedestal.first, sqrt(pedestal.second.sum2 / pedestal.second.n - (pedestal_mean * pedestal_mean)));
	}
	condIO.store("newPedestals.txt", pedestals_cond);
	condIO.store("newNoise.txt", noise_cond);

}


// ------------ method called once each job just before starting event loop  ------------
void
Pedestals::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
Pedestals::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Pedestals::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Pedestals);
