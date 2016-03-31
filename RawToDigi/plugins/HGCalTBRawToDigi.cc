#include <iostream>
#include "HGCal/RawToDigi/plugins/HGCalTBRawToDigi.h"
using namespace std;

unsigned int gray_to_binary (unsigned int gray);

HGCalTBRawToDigi::HGCalTBRawToDigi(edm::ParameterSet const& conf):
	dataTag_(conf.getParameter<edm::InputTag>("InputLabel")),
	fedId_(conf.getUntrackedParameter<int>("fedId")),
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

	std::auto_ptr<SKIROC2DigiCollection> digis(0);
	//
	const FEDRawData& fed = rawraw->FEDData(fedId_);
	if (fed.size() != 0) { /// \todo Exception if 0????

		// we can figure out the number of samples from the size of the raw data
		int nsamples = fed.size() / (sizeof(uint16_t) * SKIROC::NCHANNELS * 2); // 2 is for ADC and TDC
		digis = std::auto_ptr<SKIROC2DigiCollection>(new SKIROC2DigiCollection(nsamples));
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
                                        digis->backDataFrame().setSample(0, gray_to_binary(pdata[ptradc1] & 0xFFF),gray_to_binary( pdata[ptradc2] & 0xFFF)); 
				}
			}
                 }

	}

	// put it into the event
	e.put(digis);
}

unsigned int gray_to_binary (unsigned int gray){

  unsigned int result = gray & (1 << 11);
  result |= (gray ^ (result >> 1)) & (1 << 10);
  result |= (gray ^ (result >> 1)) & (1 << 9);
  result |= (gray ^ (result >> 1)) & (1 << 8);
  result |= (gray ^ (result >> 1)) & (1 << 7);
  result |= (gray ^ (result >> 1)) & (1 << 6);
  result |= (gray ^ (result >> 1)) & (1 << 5);
  result |= (gray ^ (result >> 1)) & (1 << 4);
  result |= (gray ^ (result >> 1)) & (1 << 3);
  result |= (gray ^ (result >> 1)) & (1 << 2);
  result |= (gray ^ (result >> 1)) & (1 << 1);
  result |= (gray ^ (result >> 1)) & (1 << 0);
  return result;
}


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(HGCalTBRawToDigi);

