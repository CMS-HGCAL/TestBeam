// -*- C++ -*-
//
// Package:    HGCal/RecHitPlotter_HighGain_Correlation_CM
// Class:      RecHitPlotter_HighGain_Correlation_CM
//
/**\class RecHitPlotter_HighGain_Correlation_CM RecHitPlotter_HighGain_Correlation_CM.cc HGCal/RecHitPlotter_HighGain_Correlation_CM/plugins/RecHitPlotter_HighGain_Correlation_CM.cc

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
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
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

static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.
bool doCommonMode_CM = 1;
double return_RMS_CM(double mean_sq, double mean)
{
	return sqrt(mean_sq - mean * mean);
}

class RecHitPlotter_HighGain_Correlation_CM : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit RecHitPlotter_HighGain_Correlation_CM(const edm::ParameterSet&);
	~RecHitPlotter_HighGain_Correlation_CM();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;

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
//TH2D *Covar_hist =  new TH2D("Covar_hist","Covar_hist",nphistar_bins-1,phistar_var,nphistar_bins-1,phistar_var);
//TH2D *Correl_hist =  new TH2D("Correl_hist","Correl_hist",nphistar_bins-1,phistar_var,nphistar_bins-1,phistar_var);

	const static int NLAYERS  = 128;
	TH2Poly *h_RecHit_layer[128];
	const static int cellx = 15;
	const static int celly = 15;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
	TH1F* Full_Cell[MAXLAYERS];
	TH1F* Half_Cell[MAXLAYERS];
	TH1F* MB_Cell[MAXLAYERS];
	TH1F* Calib_Pads[MAXLAYERS];
	TH1F* Merged_Cell[MAXLAYERS];
	TH1F  *h_digi_layer_channel[MAXSKIROCS][64][MAXLAYERS];
//        TH2F  *h_digi_layer_channel_CM[2][64];
	TH1F* Sum_Cluster_ADC;
	TH1F* AllCells_Ped;
	TH1F* AllCells_CM;
	TH2F* Noise_2D_Profile;
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
RecHitPlotter_HighGain_Correlation_CM::RecHitPlotter_HighGain_Correlation_CM(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));

//Booking 2 "hexagonal" histograms to display the sum of Rechits and the Occupancy(Hit > 5 GeV) in 1 sensor in 1 layer. To include all layers soon. Also the 1D Rechits per cell in a sensor is booked here.
	AllCells_Ped = fs->make<TH1F>("AllCells_Ped", "AllCells_Ped", 500, -250, 250);
	AllCells_CM = fs->make<TH1F>("AllCells_CM", "AllCells_CM", 500, -250, 250);
	sprintf(name, "Noise_2D_Profile_Layer");
	sprintf(title, "Noise 2D Profile Layer");
	Noise_2D_Profile = fs->make<TH2F>(name, title, 2048, 0, 2048, 8000, -4000, 4000);
	for(int ILayer = 0; ILayer < MAXLAYERS; ILayer++) {
		sprintf(name, "Full_Cell_Layer_%i", ILayer);
		sprintf(title, "Full Cell Layer %i", ILayer);
		Full_Cell[ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "Half_Cell_Layer_%i", ILayer);
		sprintf(title, "Half Cell Layer %i", ILayer);
		Half_Cell[ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "MB_Cell_Layer_%i", ILayer);
		sprintf(title, "MB Cell Layer %i", ILayer);
		MB_Cell[ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "Calib_Pads_Layer_%i", ILayer);
		sprintf(title, "Calib Pads Layer %i", ILayer);
		Calib_Pads[ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "Merged_Cell_Layer_%i", ILayer);
		sprintf(title, "Merged Cell Layer %i", ILayer);
		Merged_Cell[ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		for(int ISkiroc = 1; ISkiroc <= MAXSKIROCS; ISkiroc++) {
			for(int Channel = 0; Channel < 64; Channel++) {
				sprintf(name, "Ski_%i_Channel_%i_Layer_%i", ISkiroc, Channel, ILayer);
				sprintf(title, "Ski %i Channel %i Layer %i", ISkiroc, Channel, ILayer);
				h_digi_layer_channel[ISkiroc - 1][Channel][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
				/*
				        			 sprintf(name, "Ski_%i_Channel_%i_CM",ISkiroc,Channel);
				       			         sprintf(title, "Ski %i Channel %i CM",ISkiroc,Channel);
				                		 h_digi_layer_channel_CM[ISkiroc-1][Channel] = fs->make<TH2F>(name, title, 1000,-500., 500.,1000,-500., 500.);
				*/
			}
		}
	}

}//contructor ends here


RecHitPlotter_HighGain_Correlation_CM::~RecHitPlotter_HighGain_Correlation_CM()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RecHitPlotter_HighGain_Correlation_CM::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

	using namespace edm;

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_, Rechits);
	edm::Handle<HGCalTBRecHitCollection> Rechits1;
	event.getByToken(HGCalTBRecHitCollection_, Rechits1);

	double Average_Pedestal_Per_Event_Full[MAXLAYERS] = {0};
	int Cell_counter[MAXLAYERS] = {0};
	double Average_Pedestal_Per_Event_Half[MAXLAYERS] = {0};
	int Cell_counter_Half[MAXLAYERS] = {0};
	double Average_Pedestal_Per_Event_MB[MAXLAYERS] = {0};
	int Cell_counter_MB[MAXLAYERS] = {0};
	double Average_Pedestal_Per_Event_Calib_Pad[MAXLAYERS] = {0};
	int Cell_counter_Calib_Pad[MAXLAYERS] = {0};
	double Average_Pedestal_Per_Event_Merged_Cell[MAXLAYERS] = {0};
	int Cell_counter_Merged_Cell[MAXLAYERS] = {0};
	for(auto RecHit1 : *Rechits1) {
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit1.id()).layer(), (RecHit1.id()).sensorIU(), (RecHit1.id()).sensorIV(), (RecHit1.id()).iu(), (RecHit1.id()).iv(), sensorsize);
		uint32_t EID = essource_.emap_.detId2eid(RecHit1.id());
		HGCalTBElectronicsId eid(EID);
//                if((eid.iskiroc()%2 == 1) && (eid.ichan() == 0 || eid.ichan() == 1 )) continue;
//             double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
//             double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
//               Cell_counter++;
//               Average_Pedestal_Per_Event_Full += RecHit1.energyHigh();
		if(RecHit1.energyHigh() > 32.) continue;
		if((RecHit1.id()).cellType() == 0) {
//                       Full_Cell[(RecHit1.id()).layer() - 1]->Fill(RecHit1.energyHigh());
			Cell_counter[(RecHit1.id()).layer() - 1]++;
			Average_Pedestal_Per_Event_Full[(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		} else if((RecHit1.id()).cellType() == 2) {
//		       Half_Cell[(RecHit1.id()).layer() - 1]->Fill(RecHit1.energyHigh());
			Cell_counter_Half[(RecHit1.id()).layer() - 1]++;
			Average_Pedestal_Per_Event_Half[(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		} else if(((RecHit1.id()).cellType() == 3) && ( ((RecHit1.id()).iu() == 4 && (RecHit1.id()).iv() == 3) || ((RecHit1.id()).iu() == -7 && (RecHit1.id()).iv() == 4) || ((RecHit1.id()).iu() == 7 && (RecHit1.id()).iv() == -3) || ((RecHit1.id()).iu() == -4 && (RecHit1.id()).iv() == -3) )) {
//		     MB_Cell[(RecHit1.id()).layer() - 1]->Fill(RecHit1.energyHigh());
			Cell_counter_MB[(RecHit1.id()).layer() - 1]++;
			Average_Pedestal_Per_Event_MB[(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		} else if(((RecHit1.id()).cellType() == 3) && (((RecHit1.id()).iu() == -4 && (RecHit1.id()).iv() == 6) || ((RecHit1.id()).iu() == -2 && (RecHit1.id()).iv() == 6) || ((RecHit1.id()).iu() == 4 && (RecHit1.id()).iv() == -7) || ((RecHit1.id()).iu() == 2 && (RecHit1.id()).iv() == -6))) {
//                     Merged_Cell[(RecHit1.id()).layer() - 1]->Fill(RecHit1.energyHigh());
			Cell_counter_Merged_Cell[(RecHit1.id()).layer() - 1]++;
			Average_Pedestal_Per_Event_Merged_Cell[(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		} else if((RecHit1.id()).cellType() == 1) {
//                      Calib_Pads[(RecHit1.id()).layer() - 1]->Fill(RecHit1.energyHigh());
			Cell_counter_Calib_Pad[(RecHit1.id()).layer() - 1]++;
			Average_Pedestal_Per_Event_Calib_Pad[(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		}

	}


	for(int iii = 0; iii < MAXLAYERS; iii++) {
		if(Cell_counter[iii] != 0) Full_Cell[iii]->Fill(Average_Pedestal_Per_Event_Full[iii] / Cell_counter[iii]);
		if(Cell_counter_Half[iii] != 0) Half_Cell[iii]->Fill(Average_Pedestal_Per_Event_Half[iii] / Cell_counter_Half[iii]);
		if(Cell_counter_MB[iii] != 0) MB_Cell[iii]->Fill(Average_Pedestal_Per_Event_MB[iii] / Cell_counter_MB[iii]);
		if(Cell_counter_Merged_Cell[iii] != 0) Merged_Cell[iii]->Fill(Average_Pedestal_Per_Event_Merged_Cell[iii] / Cell_counter_Merged_Cell[iii]);
		if(Cell_counter_Calib_Pad[iii] != 0) Calib_Pads[iii]->Fill(Average_Pedestal_Per_Event_Calib_Pad[iii] / Cell_counter_Calib_Pad[iii]);
//            cout<<endl<<"iii= "<<iii<<" "<<Cell_counter[iii]<<" "<<Cell_counter_Half[iii]<<" "<<Cell_counter_MB[iii]<<" "<<Cell_counter_Calib_Pad[iii]<<" "<<Cell_counter_Merged_Cell[iii]<<endl;
	}

	for(auto RecHit : *Rechits) {
//                if(RecHit.energyHigh() > 32.) continue;
		if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize);
//                double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
//                double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
		uint32_t EID = essource_.emap_.detId2eid(RecHit.id());
		HGCalTBElectronicsId eid(EID);
//                          TF1* Fit2= (TF1*) h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()]->GetFunction("gaus");
		if(!doCommonMode_CM) {
			h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh());
			Noise_2D_Profile->Fill((64 * (eid.iskiroc() - 1) + eid.ichan()), RecHit.energyHigh());
		}
		if(doCommonMode_CM) {
			AllCells_Ped->Fill(RecHit.energyHigh());
//                             cout<<endl<<" Energy= "<<RecHit.energyHigh()<<" CM = "<<Average_Pedestal_Per_Event_Full[(RecHit.id()).layer() - 1]<<" Cells= "<<Cell_counter[(RecHit.id()).layer() -1]<<endl;
			if((RecHit.id()).cellType() == 0 || (RecHit.id()).cellType() == 4) h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_Full[(RecHit.id()).layer() - 1] / (Cell_counter[(RecHit.id()).layer() - 1])));
			else if((RecHit.id()).cellType() == 2 ) h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_Half[(RecHit.id()).layer() - 1] / (Cell_counter_Half[(RecHit.id()).layer() - 1])));
			else if((RecHit.id()).cellType() == 1) h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_Calib_Pad[(RecHit.id()).layer() - 1] / (Cell_counter_Calib_Pad[(RecHit.id()).layer() - 1])));
			else if(((RecHit.id()).cellType() == 3) && ( ((RecHit.id()).iu() == 4 && (RecHit.id()).iv() == 3) || ((RecHit.id()).iu() == -7 && (RecHit.id()).iv() == 4) || ((RecHit.id()).iu() == 7 && (RecHit.id()).iv() == -3) || ((RecHit.id()).iu() == -4 && (RecHit.id()).iv() == -3) ) ) {
				h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_MB[(RecHit.id()).layer() - 1] / (Cell_counter_MB[(RecHit.id()).layer() - 1])));
			} else if(((RecHit.id()).cellType() == 3) && (((RecHit.id()).iu() == -4 && (RecHit.id()).iv() == 6) || ((RecHit.id()).iu() == -2 && (RecHit.id()).iv() == 6) || ((RecHit.id()).iu() == 4 && (RecHit.id()).iv() == -7) || ((RecHit.id()).iu() == 2 && (RecHit.id()).iv() == -6))) {
				h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_Merged_Cell[(RecHit.id()).layer() - 1] / (Cell_counter_Merged_Cell[(RecHit.id()).layer() - 1])));
			}
			Noise_2D_Profile->Fill((64 * (eid.iskiroc() - 1) + eid.ichan()), RecHit.energyHigh() - (Average_Pedestal_Per_Event_Full[(RecHit.id()).layer() - 1] / (Cell_counter[(RecHit.id()).layer() - 1])));
		}


	}


}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
RecHitPlotter_HighGain_Correlation_CM::beginJob()
{
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(mapfile_);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void
RecHitPlotter_HighGain_Correlation_CM::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitPlotter_HighGain_Correlation_CM::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitPlotter_HighGain_Correlation_CM);
