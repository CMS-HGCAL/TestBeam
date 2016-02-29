#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"


#include "HGCal/DataFormats/interface/SKIROCParameters.h"

#include <iostream>

class HGCalTBRawToDigi : public edm::EDProducer
{
public:
	explicit HGCalTBRawToDigi(const edm::ParameterSet& ps);
	virtual void produce(edm::Event& e, const edm::EventSetup& c);
	virtual void beginJob();

private:
	edm::InputTag dataTag_;
	int fedId_;
	std::string mapfile_;
	struct {
		HGCalElectronicsMap emap_;
	} essource_;

};

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
		int skiroc = 1; // currently we don't have a way to tell these apart
		int ptr = fed.size() / sizeof(uint16_t) - 1;

		printf("Starting on SKIROC %x\n", pdata[ptr]);
		ptr--; // now we are pointing at a relatively-useless header word
		ptr--; // now we are pointing at the first TDC word

		for (int ichan = 0; ichan < SKIROC::NCHANNELS; ichan++) {
			HGCalTBElectronicsId eid(skiroc, ichan);
			if (!essource_.emap_.existsEId(eid.rawId())) {
				std::cout << "We do not have a mapping for " << eid << std::endl;
			} else {
				HGCalTBDetId did = essource_.emap_.eid2detId(eid);
				digis->addDataFrame(did);
				for (int is = 0; is < nsamples; is++) {
					int ptrtdc = ptr - ichan - 128 * is;
					int ptradc = ptr - ichan - 64 - 128 * is;
					digis->backDataFrame().setSample(is, pdata[ptradc], pdata[ptrtdc]);
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

