// -*- C++ -*-
//
// Package:    HGCal/RecHitPlotter_HighGain_New
// Class:      RecHitPlotter_HighGain_New
//
/**\class RecHitPlotter_HighGain_New RecHitPlotter_HighGain_New.cc HGCal/RecHitPlotter_HighGain_New/plugins/RecHitPlotter_HighGain_New.cc

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
#include <fstream>
#include "TH2Poly.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBTrackCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBTrack.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"
#define MAXVERTICES 6

using namespace std;

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.
double Return_RMS(double mean_sq, double mean)
{
	return sqrt(mean_sq - mean * mean);
}

edm::Service<TFileService> fs;
class RecHitPlotter_HighGain_New : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit RecHitPlotter_HighGain_New(const edm::ParameterSet&);
	~RecHitPlotter_HighGain_New();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;
	void InitTH2Poly(TH2Poly& poly, int iLayer);//Draws cells into the TH2Poly structure
	// ----------member data ---------------------------
	edm::EDGetToken HGCalTBRecHitCollection_;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
	struct {
		HGCalElectronicsMap emap_;
	} essource_;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	char name[50], title[50];
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
RecHitPlotter_HighGain_New::RecHitPlotter_HighGain_New(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
}//contructor ends here


RecHitPlotter_HighGain_New::~RecHitPlotter_HighGain_New()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

//Creating a Sensor with cells per layer////////////////////////////
void RecHitPlotter_HighGain_New::InitTH2Poly(TH2Poly& poly, int iLayer){

	double HexX[MAXVERTICES] = {0.};
	double HexY[MAXVERTICES] = {0.};

	for(int ISkiroc = (2*iLayer) - 1; ISkiroc <= (2*iLayer); ISkiroc++) {
		for(int Channel = 0; Channel < 64; Channel++) {
                	HGCalTBElectronicsId ElId(ISkiroc,Channel);
                        HGCalTBDetId DetId = essource_.emap_.eid2detId(ElId);
			CellXY = TheCell.GetCellCoordinatesForPlots(DetId.layer(), DetId.sensorIU(), DetId.sensorIV(), DetId.iu(), DetId.iv(), DetId.cellType(), sensorsize);
                        assert(CellXY.size() == 4 || CellXY.size() == 6);
                        unsigned int iVertex = 0;
                        for(it = CellXY.begin(); it != CellXY.end(); it++) {
                                HexX[iVertex] =  it->first;
                                HexY[iVertex] =  it->second;
                                ++iVertex;
                            }// loop over the vertices of a single cell
			poly.AddBin(CellXY.size(), HexX, HexY);
		     }//loop over Channels in a skiroc
	     }//loop over skirocs in a layer
      }
//////////////////////////////////////////////////////////////////////////

// ------------ method called for each event  ------------
void
RecHitPlotter_HighGain_New::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

	using namespace edm;
	using namespace std;

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_, Rechits);
//////////////////////////////Declaring the TH2Poly for each layer///////////
	TH2Poly *h_RecHit_layer[MAXLAYERS];
	int evId = event.id().event() - 1;
	int iLayer = (evId % (MAXLAYERS * EVENTSPERSPILL)) / (EVENTSPERSPILL);
	h_RecHit_layer[iLayer] = fs->make<TH2Poly>();
	sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", 0, iLayer + 1, evId);
	sprintf(title, "ADC counts in Layer%i", iLayer + 1);
	h_RecHit_layer[iLayer]->SetName(name);
	h_RecHit_layer[iLayer]->SetTitle(title);
	InitTH2Poly(*h_RecHit_layer[iLayer], iLayer + 1);
//////////////////////////////////////////////////////////////
//
	for(auto RecHit : *Rechits) {
		if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
		int n_layer = (RecHit.id()).layer();
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), (RecHit.id()).cellType(), sensorsize);
		double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
		double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
		h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, RecHit.energyHigh());//Filling a TH2Poly
	}

}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
RecHitPlotter_HighGain_New::beginJob()
{
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(mapfile_);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void
RecHitPlotter_HighGain_New::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitPlotter_HighGain_New::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitPlotter_HighGain_New);
