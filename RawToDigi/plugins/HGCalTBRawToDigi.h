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

/**
 * \class HGCal/RawToDigi/plugins/HGCalTBRawToDigi.h HGCalTBRawToDigi.h HGCalTBRawToDigi
 * \brief Produces a digi collection starting from FEDRawData
 */


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
	int ptradc1, ptradc2;// location HG,LG SKI,Channel in the pushed vector of 16 bit words.
};
