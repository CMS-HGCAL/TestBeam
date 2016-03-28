#include <iostream>
#include "HGCal/RawToDigi/plugins/HGCalTBRawToDigi.h"
using namespace std;

HGCalTBRawToDigi::HGCalTBRawToDigi(edm::ParameterSet const& conf):
	dataTag_(conf.getParameter<edm::InputTag>("InputLabel")),
	fedIds_(conf.getUntrackedParameter<std::vector<int> >("fedIds")),
	mapfile_(conf.getUntrackedParameter<std::string>("electronicsMap"))
{
	produces<SKIROC2DigiCollection>();
	consumes<FEDRawDataCollection>(dataTag_);
}

void HGCalTBRawToDigi::beginJob()
{
	///\todo this should be done by an ESProducer
	HGCalCondObjectTextIO io(0); // don't need a numbering scheme for this

	edm::FileInPath fip(mapfile_);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	}
}

void HGCalTBRawToDigi::produce(edm::Event& e, const edm::EventSetup& c)
{
	// Step A: Get Inputs
	edm::Handle<FEDRawDataCollection> rawraw;
	e.getByLabel(dataTag_, rawraw);

	std::auto_ptr<SKIROC2DigiCollection> digis = std::auto_ptr<SKIROC2DigiCollection>(new SKIROC2DigiCollection(SKIROC::MAXSAMPLES)); // 2 because TDC and ADC
	//

	for(auto fedId_ : fedIds_){
	  int skiroc = fedId_;
	  std::cout << skiroc << std::endl;
	  const FEDRawData& fed = rawraw->FEDData(fedId_);
	  if(fed.size() == 0) continue; // empty FEDs are allowed: not in the readout for example

		// we can figure out the number of samples from the size of the raw data
		//int nsamples = fed.size() / (sizeof(uint16_t) * SKIROC::NCHANNELS * 2); // 2 is for ADC and TDC
		//assert(nsamples>0);
		
		const uint16_t* pdata = (const uint16_t*)(fed.data());

		// we start from the back...
		int ptr = fed.size() / sizeof(uint16_t) - 1;
/*
		printf("Starting on SKIROC %x\n", pdata[ptr]);
		ptr--; // now we are pointing at a relatively-useless header word
		ptr--; // now we are pointing at the first TDC word
*/
                for (int ski = 2; ski >= 1; ski--) {
			for (int ichan = 0; ichan < SKIROC::NCHANNELS; ichan++) {
				HGCalTBElectronicsId eid(ski, ichan);
				if (!essource_.emap_.existsEId(eid.rawId())) {
					std::cout << "We do not have a mapping for " << eid << std::endl;
				} else {
					HGCalTBDetId did = essource_.emap_.eid2detId(eid);
					digis->addDataFrame(did);
					int ptradc1 = ptr - ichan - 128 * (2 - ski);
					int ptradc2 = ptr - ichan - 64 - 128 * (2 - ski);
//					digis->backDataFrame().setSample(0, pdata[ptradc1], pdata[ptradc2]);// Only one sample hardcoded as 0
                                        digis->backDataFrame().setSample(0, pdata[ptradc2], pdata[ptradc1]);  
				}
			}
                 }

	}


	// put it into the event
	e.put(digis);
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(HGCalTBRawToDigi);

