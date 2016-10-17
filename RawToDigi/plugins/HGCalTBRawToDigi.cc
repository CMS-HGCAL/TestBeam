#include "HGCal/RawToDigi/plugins/HGCalTBRawToDigi.h"

// provides the maximum number of layers and number of skirocs per board
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

// provide the number of channels in a skiroc and number of samples
#include "HGCal/DataFormats/interface/SKIROCParameters.h"

//#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif

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
		size_t nSkirocs = fed.size() / (sizeof(uint16_t) * SKIROC::NCHANNELS * SKIROC::MAXSAMPLES);
		size_t nBoards = nSkirocs / MAXSKIROCS_PER_BOARD;
		digis = std::auto_ptr<SKIROC2DigiCollection>(new SKIROC2DigiCollection(nSkirocs * SKIROC::NCHANNELS * SKIROC::MAXSAMPLES));
		const uint16_t* pdata = (const uint16_t*)(fed.data());

		size_t ski = 0; // the skirocs have an absolute numbering, start counting from the first board till the last
		for (unsigned int i_board = 0 ; i_board < nBoards; ++i_board) {
			for(size_t i_skiroc = 0; i_skiroc < MAXSKIROCS_PER_BOARD; ++i_skiroc) {
				for (int ichan = 0; ichan < SKIROC::NCHANNELS; ichan++) {
					HGCalTBElectronicsId eid(ski, ichan);
					if (essource_.emap_.existsEId(eid.rawId())) {
						HGCalTBDetId did = essource_.emap_.eid2detId(eid);
						digis->addDataFrame(did);
#ifdef DEBUG
						if(i_board == 0) std::cout << (*pdata & 0xFFF) << "\t" << (*(pdata+1) & 0xFFF) << "\t" << (*(pdata+2) & 0xFFF) << std::endl;
#endif
						digis->backDataFrame().setSample(0, gray_to_binary(*(pdata) & 0xFFF), gray_to_binary( *(pdata+1) & 0xFFF), 0);
						pdata++;pdata++;
					}
				}
				++ski; //increment the absolute ID of the skiroc
			}
		}


	}// fed size > 0

	// put it into the event
	e.put(digis);
}

unsigned int gray_to_binary (unsigned int gray)
{

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

