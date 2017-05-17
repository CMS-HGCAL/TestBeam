// -*- C++ -*-
//
// Package:    HGCal/DigiPlotter
// Class:      DigiPlotter
//
/**\class DigiPlotter DigiPlotter.cc HGCal/DigiPlotter/plugins/DigiPlotter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Mon, 15 Feb 2016 09:47:43 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TH1F.h"
#include "TProfile.h"
#include <fstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

using namespace std;


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DigiPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit DigiPlotter(const edm::ParameterSet&);
	~DigiPlotter();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;
	// ----------member data ---------------------------
	bool DEBUG = 0;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt";
	struct {
		HGCalElectronicsMap emap_;
	} essource_;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	const static int NSAMPLES = 2;
	char name[50], title[50];
};

DigiPlotter::DigiPlotter(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	consumesMany<HGCalTBSkiroc2CMSCollection>();

}//contructor ends here


DigiPlotter::~DigiPlotter()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DigiPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	using namespace edm;
	std::vector<edm::Handle<HGCalTBSkiroc2CMSCollection>> ski;
	event.getManyByType(ski);



}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
DigiPlotter::beginJob()
{
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(mapfile_);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void
DigiPlotter::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DigiPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	descriptions.add("hgcaltbdigisplotter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiPlotter);
