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
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"

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
	bool DEBUG = 1;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	const static int NSAMPLES = 1;
	const static int NLAYERS  = 4;
	TH2Poly *h_digi_layer[NSAMPLES][NLAYERS];
        TH1F    *h_digi_layer_summed[NSAMPLES][NLAYERS];
	const static int cellx = 15;
	const static int celly = 15;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
	TH1F  *h_digi_layer_cell[NSAMPLES][NLAYERS][cellx][celly];
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
DigiPlotter::DigiPlotter(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	consumesMany<SKIROC2DigiCollection>();
	const int HalfHexVertices = 4;
	double HalfHexX[HalfHexVertices] = {0.};
	double HalfHexY[HalfHexVertices] = {0.};
	const int FullHexVertices = 6;
	double FullHexX[FullHexVertices] = {0.};
	double FullHexY[FullHexVertices] = {0.};
	int iii = 0;
	for(int nsample = 0; nsample < NSAMPLES; nsample++) {
		for(int nlayers = 0; nlayers < NLAYERS; nlayers++) {
//Booking a "hexagonal" histograms to display the sum of Digis for NSAMPLES, in 1 SKIROC in 1 layer. To include all layers soon. Also the 1D Digis per cell in a sensor is booked here for NSAMPLES.
			sprintf(name, "FullLayer_Sample%i_Layer%i", nsample, nlayers + 1);
			sprintf(title, "Sum of adc counts per cell for Sample%i Layer%i", nsample, nlayers + 1);
			h_digi_layer[nsample][nlayers] = fs->make<TH2Poly>();
			h_digi_layer[nsample][nlayers]->SetName(name);
			h_digi_layer[nsample][nlayers]->SetTitle(title);
			sprintf(name, "FullLayer_Sample%i_Layer%i_summed", nsample, nlayers + 1);
                        sprintf(title, "Sum of adc counts for all cells in Sample%i Layer%i", nsample, nlayers + 1);
                        h_digi_layer_summed[nsample][nlayers] = fs->make<TH1F>(name, title, 4096, 0., 4095.);
                        h_digi_layer_summed[nsample][nlayers]->GetXaxis()->SetTitle("Digis[adc counts]");
			for(int iv = -7; iv < 8; iv++) {
				for(int iu = -7; iu < 8; iu++) {
					if(!IsCellValid.iu_iv_valid(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
//Some thought needs to be put in about the binning and limits of this 1D histogram, probably different for beam type Fermilab and cern.
					sprintf(name, "Cell_u_%i_v_%i_Sample%i_Layer%i", iu, iv, nsample, nlayers + 1);
					sprintf(title, "Digis for Cell_u_%i_v_%i Sample%i Layer%i", iu, iv, nsample, nlayers + 1);
					h_digi_layer_cell[nsample][nlayers][iu + 7][iv + 7] = fs->make<TH1F>(name, title, 4096, 0., 4095.);
					h_digi_layer_cell[nsample][nlayers][iu + 7][iv + 7]->GetXaxis()->SetTitle("Digis[adc counts]");
					CellXY = TheCell.GetCellCoordinatesForPlots(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
					int NumberOfCellVertices = CellXY.size();
					iii = 0;
					if(NumberOfCellVertices == 4) {
						for(it = CellXY.begin(); it != CellXY.end(); it++) {
							HalfHexX[iii] =  it->first;
							HalfHexY[iii++] =  it->second;
						}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
						h_digi_layer[nsample][nlayers]->AddBin(NumberOfCellVertices, HalfHexX, HalfHexY);
					} else if(NumberOfCellVertices == 6) {
						iii = 0;
						for(it = CellXY.begin(); it != CellXY.end(); it++) {
							FullHexX[iii] =  it->first;
							FullHexY[iii++] =  it->second;
						}
						h_digi_layer[nsample][nlayers]->AddBin(NumberOfCellVertices, FullHexX, FullHexY);
					}

				}//loop over iu
			}//loop over iv
		}//loop over nlayers
	}//loop over nsamples
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
	using namespace std;
	std::vector<edm::Handle<SKIROC2DigiCollection> > ski;
	event.getManyByType(ski);
	if(!ski.empty()) {

		std::vector<edm::Handle<SKIROC2DigiCollection> >::iterator i;
		for(i = ski.begin(); i != ski.end(); i++) {
			const SKIROC2DigiCollection& Coll = *(*i);
			cout << "SKIROC2 Digis: " << i->provenance()->branchName() << endl;
			for(SKIROC2DigiCollection::const_iterator j = Coll.begin(); j != Coll.end(); j++) {
				const SKIROC2DataFrame& SKI = *j ;
				int n_layer = (SKI.detid()).layer();
				int n_sensor_IU = (SKI.detid()).sensorIU();
				int n_sensor_IV = (SKI.detid()).sensorIV();
				int n_cell_iu = (SKI.detid()).iu();
				int n_cell_iv = (SKI.detid()).iv();
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
				if(!IsCellValid.iu_iv_valid(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize))  continue;
				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize);
				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);
				for(int nsample = 0; nsample < SKI.samples(); nsample++) {
					h_digi_layer[nsample][n_layer - 1]->Fill(iux , iyy, SKI[nsample].adc());
					h_digi_layer_summed[nsample][n_layer - 1]->Fill(SKI[nsample].adc());
					h_digi_layer_cell[nsample][n_layer - 1][7 + n_cell_iu][7 + n_cell_iv]->Fill(SKI[nsample].adc());
				}
			}
		}
	} else {
		edm::LogWarning("DQM") << "No SKIROC2 Digis";
	}

}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
DigiPlotter::beginJob()
{

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
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiPlotter);
