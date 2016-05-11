#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "HGCal/DataFormats/interface/HGCalTBTrackCollection.h"

/**
 * \class HGCal/RawToDigi/plugins/FNALTelescopeRawToTracks.h FNALTelescopeRawToTracks.h FNALTelescopeRawToTracks
 * \brief Produces a track collection starting from FEDRawData
 *
 * Unpacking done according to what implemented in
 * HGCalTBTextSource::parseAddTelescopeWords
 */


class FNALTelescopeRawToTracks : public edm::EDProducer
{
public:
	explicit FNALTelescopeRawToTracks(const edm::ParameterSet& ps);
	virtual void produce(edm::Event& e, const edm::EventSetup& c);
	virtual void beginJob();
	void fillDescription(edm::ParameterSetDescription &desc);

private:
	edm::InputTag dataTag_;
	std::vector<int> fedIds_;

};

FNALTelescopeRawToTracks::FNALTelescopeRawToTracks(edm::ParameterSet const& conf):
	dataTag_(conf.getParameter<edm::InputTag>("InputLabel")),
	fedIds_(conf.getUntrackedParameter<std::vector<int> >("fedIds"))
{
	produces<HGCalTBTrackCollection>();
	consumes<FEDRawDataCollection>(dataTag_);
}

void FNALTelescopeRawToTracks::beginJob()
{
}

void FNALTelescopeRawToTracks::produce(edm::Event& e, const edm::EventSetup& c)
{
	edm::Handle<FEDRawDataCollection> rawraw;
	e.getByLabel(dataTag_, rawraw);

	std::auto_ptr<HGCalTBTrackCollection> tracks = std::auto_ptr<HGCalTBTrackCollection>(new HGCalTBTrackCollection);

	for(auto fedId_ : fedIds_) {

		const FEDRawData& fed = rawraw->FEDData(fedId_);
		if(fed.size() == 0) continue; // empty FEDs are allowed: for example when one FED is not in the readout for a particular run

		const float* pdata = (const float*)(fed.data());

		unsigned int size = fed.size() / HGCalTBTrack::getSizeof(); // returns the number of HGCalTBTrack
		///\todo missing check if the size is correct
		for(unsigned int i = 0; i < size; ++i, pdata += HGCalTBTrack::getSizeof()) {
			HGCalTBTrack track(pdata);
			tracks->push_back(track);
		}
	}




	// put it into the event
	e.put(tracks);
}

void FNALTelescopeRawToTracks::fillDescription(edm::ParameterSetDescription &desc)
{
	desc.setComment("TEST");
	desc.add<edm::InputTag>("InputLabel");
	desc.addUntracked<std::vector<int> >("fedIds", std::vector<int>(99));
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FNALTelescopeRawToTracks);

