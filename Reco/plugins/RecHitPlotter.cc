// -*- C++ -*-
//
// Package:    HGCal/RecHitPlotter
// Class:      RecHitPlotter
//
/**\class RecHitPlotter RecHitPlotter.cc HGCal/RecHitPlotter/plugins/RecHitPlotter.cc

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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.

class RecHitPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit RecHitPlotter(const edm::ParameterSet&);
	~RecHitPlotter();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;

	// ----------member data ---------------------------
	edm::EDGetToken HGCalTBRecHitCollection_;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	const static int NLAYERS  = 2;
	TH2Poly *h_RecHit_layer[NLAYERS];
	TH1F    *h_RecHit_layer_summed[NLAYERS];
	TH2Poly *h_RecHit_layer_Occupancy[NLAYERS];
	const static int cellx = 15;
	const static int celly = 15;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
	TH1F  *h_RecHit_layer_cell[NLAYERS][cellx][celly];
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
RecHitPlotter::RecHitPlotter(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));

//Booking 2 "hexagonal" histograms to display the sum of Rechits and the Occupancy(Hit > 5 GeV) in 1 sensor in 1 layer. To include all layers soon. Also the 1D Rechits per cell in a sensor is booked here.
	const int HalfHexVertices = 4;
	double HalfHexX[HalfHexVertices] = {0.};
	double HalfHexY[HalfHexVertices] = {0.};
	const int FullHexVertices = 6;
	double FullHexX[FullHexVertices] = {0.};
	double FullHexY[FullHexVertices] = {0.};
	int iii = 0;
	for(int nlayers = 0; nlayers < NLAYERS; nlayers++) {
		sprintf(name, "FullLayer_RecHits_Layer%i", nlayers + 1);
		sprintf(title, "Sum of RecHits Layer%i", nlayers + 1);
		h_RecHit_layer[nlayers] = fs->make<TH2Poly>();
		h_RecHit_layer[nlayers]->SetName(name);
		h_RecHit_layer[nlayers]->SetTitle(title);
		sprintf(name, "FullLayer_RecHits_Layer%i_Summed", nlayers + 1);
		sprintf(title, "Sum of RecHits Layer%i Summed over the cells", nlayers + 1);
		h_RecHit_layer_summed[nlayers] = fs->make<TH1F>(name, title, 200, -100., 100.);
		h_RecHit_layer_summed[nlayers]->GetXaxis()->SetTitle("RecHits[GeV]");
		sprintf(name, "FullLayer_Occupancy_Layer%i", nlayers + 1);
		sprintf(title, "Sum of Occupancy Layer%i", nlayers + 1);
		h_RecHit_layer_Occupancy[nlayers] = fs->make<TH2Poly>();
		h_RecHit_layer_Occupancy[nlayers]->SetName(name);
		h_RecHit_layer_Occupancy[nlayers]->SetTitle(title);
		for(int iv = -7; iv < 8; iv++) {
			for(int iu = -7; iu < 8; iu++) {
				if(!IsCellValid.iu_iv_valid(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
//Some thought needs to be put in about the binning and limits of this 1D histogram, probably different for beam type Fermilab and cern.
				sprintf(name, "Cell_RecHits_u_%i_v_%i_Layer%i", iu, iv, nlayers + 1);
				sprintf(title, "RecHits for Cell_u_%i_v_%i Layer%i", iu, iv, nlayers + 1);
				h_RecHit_layer_cell[nlayers][iu + 7][iv + 7] = fs->make<TH1F>(name, title, 200, -100., 100.); // need to finalize binning
				h_RecHit_layer_cell[nlayers][iu + 7][iv + 7]->GetXaxis()->SetTitle("RecHits[GeV]");
				CellXY = TheCell.GetCellCoordinatesForPlots(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
				int NumberOfCellVertices = CellXY.size();
				iii = 0;
				if(NumberOfCellVertices == 4) {
					for(it = CellXY.begin(); it != CellXY.end(); it++) {
						HalfHexX[iii] =  it->first;
						HalfHexY[iii++] =  it->second;
					}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
					h_RecHit_layer[nlayers]->AddBin(NumberOfCellVertices, HalfHexX, HalfHexY);
					h_RecHit_layer_Occupancy[nlayers]->AddBin(NumberOfCellVertices, HalfHexX, HalfHexY);
				} else if(NumberOfCellVertices == 6) {
					iii = 0;
					for(it = CellXY.begin(); it != CellXY.end(); it++) {
						FullHexX[iii] =  it->first;
						FullHexY[iii++] =  it->second;
					}
					h_RecHit_layer[nlayers]->AddBin(NumberOfCellVertices, FullHexX, FullHexY);
					h_RecHit_layer_Occupancy[nlayers]->AddBin(NumberOfCellVertices, FullHexX, FullHexY);
				}

			}
		}
	}//loop over layers end here


}//contructor ends here


RecHitPlotter::~RecHitPlotter()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RecHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

	using namespace edm;
	using namespace std;

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_, Rechits);
	edm::Handle<HGCalTBRecHitCollection> Rechits1;
	event.getByToken(HGCalTBRecHitCollection_, Rechits1);

	double Average_Pedestal_Per_Event1_Full = 0;
	int Cell_counter1_Full = 0;
	for(auto RecHit1 : *Rechits1) {
		Average_Pedestal_Per_Event1_Full += RecHit1.energyHigh();
		Cell_counter1_Full++;
	}

	for(auto RecHit : *Rechits) {
		if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
		int n_layer = (RecHit.id()).layer();
		int n_cell_type = (RecHit.id()).cellType();
		double n_cell_area = IsCellValid.Cell_Area(n_cell_type);
//We now obtain the cartesian coordinates of the cell corresponding to an iu,iv. This may either be a full hex, a half hex or an invalid cell. If a cell is invalid based on the iu,iv index -123456 is returned for its x,y vertices
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize);
		double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
		double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
		h_RecHit_layer[n_layer - 1]->Fill(iux , iyy, RecHit.energyHigh());
		h_RecHit_layer_summed[n_layer - 1]->Fill(RecHit.energyHigh());
//The energyHigh threshold for the occupancy has been hardcoded here. Need to decide what a good choice is. Maybe dynamic per cell depending on the pedestal
		if(RecHit.energyHigh() > 5) h_RecHit_layer_Occupancy[n_layer - 1]->Fill(iux , iyy, 1. / n_cell_area);
// There will be several array indices iu, iv that wont be filled due to it being invalid. Can think of alternate array filling.
		h_RecHit_layer_cell[n_layer - 1][7 + (RecHit.id()).iu()][7 + (RecHit.id()).iv()]->Fill(RecHit.energyHigh() - Average_Pedestal_Per_Event1_Full / Cell_counter1_Full);
	}


}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
RecHitPlotter::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
RecHitPlotter::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitPlotter);
