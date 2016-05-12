// -*- C++ -*-
//
// Package:
// Class:      TestRecHit
//
/**\class TestRecHit TestRecHit.cc
   HGCal/TBStandaloneSimulatpr/plugins/TestRecHit.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Harrison B. Prosper
//         Created:  Sat, 09 Apr 2016 13:37:51 GMT


// system include files
#include <memory>
#include <iostream>
#include "TH1F.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class TestRecHit : public edm::EDAnalyzer
{
public:
	explicit TestRecHit(const edm::ParameterSet&);
	~TestRecHit();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob();
	void analyze(const edm::Event& , const edm::EventSetup&);
	virtual void endJob();

	// ----------member data ---------------------------
	int sensor_size;
	edm::EDGetToken token;
	HGCalTBTopology check;
};

TestRecHit::TestRecHit(const edm::ParameterSet& config)
	: sensor_size(128),
	  token(consumes<HGCalTBRecHitCollection>
	        (config.getParameter<edm::InputTag>("HGCALTBRECHITS"))),
	  check(HGCalTBTopology())
{
}


TestRecHit::~TestRecHit()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TestRecHit::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	using namespace edm;
	using namespace std;

	Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(token, Rechits);

	for(auto RecHit : *Rechits) {
		int layer    = RecHit.id().layer();
		int sensor_u = RecHit.id().sensorIU();
		int sensor_v = RecHit.id().sensorIV();
		int u        = RecHit.id().iu();
		int v        = RecHit.id().iv();
		bool good = check.iu_iv_valid(layer,
		                              sensor_u, sensor_v,
		                              u, v,
		                              sensor_size);
		if ( !good ) {
			char record[1024];
			sprintf(record,
			        "bad cell:\n"
			        "\tlayer=%4d, sensor_u=%4d, sensor_v=%4d, u=%4d, v=%4d",
			        layer, sensor_u, sensor_v, u, v);
			cout << record << endl;
			continue;
		}
	}
}


void
TestRecHit::beginJob()
{
}

void
TestRecHit::endJob()
{
}

void
TestRecHit::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(TestRecHit);
