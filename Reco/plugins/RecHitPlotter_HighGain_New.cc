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
bool DoCommonMode = 1;
int Event_multiple = 1;

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
	void InitTH2Poly(TH2Poly& poly);
	ofstream ofs;
	// ----------member data ---------------------------
	edm::EDGetToken HGCalTBRecHitCollection_;
	edm::EDGetToken HGCalTBTrackCollection_;
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
	const static int NLAYERS  = 1;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
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
//Booking 2 "hexagonal" histograms to display the sum of Rechits and the Occupancy(Hit > 5 GeV) in 1 sensor in 1 layer. To include all layers soon. Also the 1D Rechits per cell in a sensor is booked here.
}//contructor ends here


RecHitPlotter_HighGain_New::~RecHitPlotter_HighGain_New()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void RecHitPlotter_HighGain_New::InitTH2Poly(TH2Poly& poly)
{
	double HexX[MAXVERTICES] = {0.};
	double HexY[MAXVERTICES] = {0.};

	for(int iv = -7; iv < 8; iv++) {
		for(int iu = -7; iu < 8; iu++) {
			if(!IsCellValid.iu_iv_valid(NLAYERS, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
			CellXY = TheCell.GetCellCoordinatesForPlots(NLAYERS, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
			assert(CellXY.size() == 4 || CellXY.size() == 6);
			unsigned int iVertex = 0;
			for(it = CellXY.begin(); it != CellXY.end(); it++) {
				HexX[iVertex] =  it->first;
				HexY[iVertex] =  it->second;
				++iVertex;
			}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
			poly.AddBin(CellXY.size(), HexX, HexY);
		}//loop over iu
	}//loop over iv
}
//

// ------------ method called for each event  ------------
void
RecHitPlotter_HighGain_New::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

	using namespace edm;
	using namespace std;

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_, Rechits);
	edm::Handle<HGCalTBRecHitCollection> Rechits1;
	event.getByToken(HGCalTBRecHitCollection_, Rechits1);

	double Average_Pedestal_Per_Event_Full = 0, Average_Pedestal_Per_Event_Half = 0, Average_Pedestal_Per_Event_MB = 0, Average_Pedestal_Per_Event_Merged = 0;
	int Cell_counter_Full = 0, Cell_counter_Half = 0, Cell_counter_MB = 0, Cell_counter_Merged = 0;

	/*
		for(auto Track : *Tracks) {
			cout << endl << " Event/Trigger = " << event.id().event() << " Track intercept X= " << Track.vertex().X() << " Track intercept Y= " << Track.vertex().Y() << " Slope X= " << Track.momentum().X() << " Slope Y= " << Track.momentum().Y() << " Track hit X= " << Track.pointAt(ExtrapolateZ).X() << " Track hit Y= " << Track.pointAt(ExtrapolateZ).Y() << endl;
		}
	*/

	for(auto RecHit1 : *Rechits1) {
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit1.id()).layer(), (RecHit1.id()).sensorIU(), (RecHit1.id()).sensorIV(), (RecHit1.id()).iu(), (RecHit1.id()).iv(), sensorsize);

		if(RecHit1.energyHigh() > 30.) continue;

		if(((RecHit1.id()).cellType() == 0) || ((RecHit1.id()).cellType() == 4)) {

				Cell_counter_Full++;
				Average_Pedestal_Per_Event_Full += RecHit1.energyHigh();

		}

		if((RecHit1.id()).cellType() == 5) {

                                Cell_counter_Merged++;
                                Average_Pedestal_Per_Event_Merged += RecHit1.energyHigh();
                        
                }

		if((RecHit1.id()).cellType() == 3) {

                                Cell_counter_MB++;
                                Average_Pedestal_Per_Event_MB += RecHit1.energyHigh();

                }

		if((RecHit1.id()).cellType() == 2) {

				Cell_counter_Half++;
				Average_Pedestal_Per_Event_Half += RecHit1.energyHigh();
			
		}

	}
	TH2Poly *h_RecHit_layer[MAXLAYERS];
	int evId = event.id().event();
	int spillId = event.luminosityBlock();
	if(((evId - 1)%Event_multiple) == 0){
		for(int iLayer = 1; iLayer<= MAXLAYERS; iLayer++){
			h_RecHit_layer[iLayer - 1] = fs->make<TH2Poly>();
			sprintf(name, "FullLayer_ADC%i_Layer%i_Spill%i_Event%i", 0, iLayer, spillId, evId);
			sprintf(title, "ADC counts in Layer%i Spill%i Event%i", iLayer, spillId, evId);
			h_RecHit_layer[iLayer -1]->SetName(name);
			h_RecHit_layer[iLayer -1]->SetTitle(title);
			InitTH2Poly(*h_RecHit_layer[iLayer -1]);
		}		
	}

	

	for(auto RecHit : *Rechits) {
		if(((evId - 1)%Event_multiple) != 0) continue;
		
		if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
		int n_layer = (RecHit.id()).layer();
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize);
		double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
		double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
		uint32_t EID = essource_.emap_.detId2eid(RecHit.id());
		HGCalTBElectronicsId eid(EID);
		if(!DoCommonMode) {
			h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, RecHit.energy());
		}
		else{//If Common mode subtraction is enabled
			if(((RecHit.id()).cellType() == 0 ) || ((RecHit.id()).cellType() == 4) ) {
				h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, (RecHit.energy() - (Average_Pedestal_Per_Event_Full / (Cell_counter_Full))) );
			}
			if((RecHit.id()).cellType() == 1 ) {
				h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, RecHit.energy());
			}

			if((RecHit.id()).cellType() == 2 ) {
				h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, RecHit.energy() - (Average_Pedestal_Per_Event_Half / Cell_counter_Half) );
			}
			if((RecHit.id()).cellType() == 3 ) {
                                h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, RecHit.energy() - (Average_Pedestal_Per_Event_MB/Cell_counter_MB) );
                        }
			if((RecHit.id()).cellType() == 5 ) {
                                h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, RecHit.energy() - (Average_Pedestal_Per_Event_Merged/Cell_counter_Merged) );
                        }

		}//else common mode condition
	}//Rechits loop ends here


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
	ofs.open("ShowerMax_X_Y_Ped_20532.txt", std::ofstream::out);
}

// ------------ method called once each job just after ending the event loop  ------------
void
RecHitPlotter_HighGain_New::endJob()
{

	ofs.close();

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
