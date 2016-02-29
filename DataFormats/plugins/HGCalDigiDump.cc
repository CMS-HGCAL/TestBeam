#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"

#include <iostream>

using namespace std;

/** \class HGCalDigiDump

\author J. Mans - Minnesota
*/
class HGCalDigiDump : public edm::EDAnalyzer
{
public:
	explicit HGCalDigiDump(edm::ParameterSet const& conf);
	virtual void analyze(edm::Event const& e, edm::EventSetup const& c);
};


HGCalDigiDump::HGCalDigiDump(edm::ParameterSet const& conf)
{
	consumesMany<SKIROC2DigiCollection>();
}
void HGCalDigiDump::analyze(edm::Event const& e, edm::EventSetup const& c)
{

	std::vector<edm::Handle<SKIROC2DigiCollection> > ski;

	try {
		e.getManyByType(ski);
		std::vector<edm::Handle<SKIROC2DigiCollection> >::iterator i;
		for (i = ski.begin(); i != ski.end(); i++) {
			const SKIROC2DigiCollection& c = *(*i);

			cout << "SKIROC2 Digis: " << i->provenance()->branchName() << endl;

			for (SKIROC2DigiCollection::const_iterator j = c.begin(); j != c.end(); j++)
				cout << "  " << *j << std::endl;
		}
	} catch (...) {
		cout << "No SKIROC2 Digis." << endl;
	}
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(HGCalDigiDump);
