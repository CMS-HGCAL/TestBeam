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
bool DoCommonMode = 0;
bool PED = 0;
bool UP = 1;

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
	ofstream ofs;
	// ----------member data ---------------------------
	edm::EDGetToken HGCalTBRecHitCollection_;
	edm::EDGetToken HGCalTBTrackCollection_;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	std::string mapfile_ = "HGCal/CondObjects/data/map_FNAL_2.txt";
	struct {
		HGCalElectronicsMap emap_;
	} essource_;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	double Mean_SUM[128] = {0.};
	double Mean_SQ_SUM[128] = {0.};
	double Mean_PROD_SUM[128][128] = { {0.} };
	double Diff_IJ_SUM[128][128] = { {0.} };
	double Mean_i = 0.;
	double Mean_j = 0.;
	double RMS_i = 0.;
	double RMS_j = 0.;
	double Correlation = 0.;
	double Covariance = 0.;
	unsigned long  Time_tmp = 0;
	unsigned long  Time_int = 0;
	TH1F* Time_hist;
	TH2F *Covar_hist;
	TH2F *Correl_hist;
	TH2F *DiffIJ_hist;
//TH2D *Covar_hist =  new TH2D("Covar_hist","Covar_hist",nphistar_bins-1,phistar_var,nphistar_bins-1,phistar_var);
//TH2D *Correl_hist =  new TH2D("Correl_hist","Correl_hist",nphistar_bins-1,phistar_var,nphistar_bins-1,phistar_var);

	const static int NLAYERS  = 1;
	TH2Poly *h_RecHit_layer[NLAYERS][1000];
	const static int cellx = 15;
	const static int celly = 15;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
	double ADC_Chan[128][9000];
	TH1F    *h_RecHit_layer_summed[NLAYERS];
	TH1F  *h_RecHit_layer_cell[NLAYERS][cellx][celly];
	TH1F* Sum_Cluster_ADC;
	TH1F* Sum_Cluster_Max;
	TH1F* AllCells_Ped;
	TH1F* AllCells_CM;
	TProfile* SKI1_Ped_Event;
	TProfile* SKI2_Ped_Event;
	TH1F* CG_X;
	TH1F* CG_Y;
	char name[50], title[50];
	double ExtrapolateZ = 20.;// distance along Z the tracks are extrapolated by(in cm) to hit the first layer.
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
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	HGCalTBTrackCollection_  = consumes<HGCalTBTrackCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBTRACKS"));
//Booking 2 "hexagonal" histograms to display the sum of Rechits and the Occupancy(Hit > 5 GeV) in 1 sensor in 1 layer. To include all layers soon. Also the 1D Rechits per cell in a sensor is booked here.
	const int HalfHexVertices = 4;
	double HalfHexX[HalfHexVertices] = {0.};
	double HalfHexY[HalfHexVertices] = {0.};
	const int FullHexVertices = 6;
	double FullHexX[FullHexVertices] = {0.};
	double FullHexY[FullHexVertices] = {0.};
	int iii = 0;
	Time_hist = fs->make<TH1F>("Time_hist", "Time_hist", 500, 1, 100000.);
	Covar_hist  = fs->make<TH2F>("Covar_hist", "Covar_hist", 128, 0, 128, 128, 0, 128);
	Correl_hist = fs->make<TH2F>("Correl_hist", "Correl_hist", 128, 0, 128, 128, 0, 128);
	DiffIJ_hist = fs->make<TH2F>("DiffIJ_hist", "DiffIJ_hist", 128, 0, 128, 128, 0, 128);
	CG_X = fs->make<TH1F>("CG_X", "CG X[cm]", 100, -5, 5);
	CG_Y = fs->make<TH1F>("CG_Y", "CG Y[cm]", 100, -5, 5);
	AllCells_Ped = fs->make<TH1F>("AllCells_Ped", "AllCells_Ped", 500, -250, 250);
	AllCells_CM = fs->make<TH1F>("AllCells_CM", "AllCells_CM", 500, -250, 250);
	SKI1_Ped_Event = fs->make<TProfile>("SKI1_Ped_Event", "Profile per event of pedestal for SKI1", 5000, 1, 5000, 0, 400);
	SKI2_Ped_Event = fs->make<TProfile>("SKI2_Ped_Event", "Profile per event of pedestal for SKI2", 5000, 1, 5000, 0, 400);
	Sum_Cluster_ADC = fs->make<TH1F>("Sum_Cluster_ADC", "Sum_Cluster_ADC",  1000, -4000., 4000.);
	Sum_Cluster_Max = fs->make<TH1F>("Sum_Cluster_Max", "Sum_Cluster_Max",  1000, -4000., 4000.);
	for(int nlayers = 0; nlayers < NLAYERS; nlayers++) {
		sprintf(name, "FullLayer_RecHits_Layer%i_Summed", nlayers + 1);
		sprintf(title, "Sum of RecHits Layer%i Summed over the cells", nlayers + 1);
		h_RecHit_layer_summed[nlayers] = fs->make<TH1F>(name, title, 4000, -2000., 2000.);
		h_RecHit_layer_summed[nlayers]->GetXaxis()->SetTitle("RecHits[GeV]");
		for(int eee = 0; eee < 1000; eee++) {
			sprintf(name, "FullLayer_RecHits_Layer%i_Event%i", nlayers + 1, eee);
			sprintf(title, "RecHits Layer%i Event%i", nlayers + 1, eee);
			h_RecHit_layer[nlayers][eee] = fs->make<TH2Poly>();
			h_RecHit_layer[nlayers][eee]->SetName(name);
			h_RecHit_layer[nlayers][eee]->SetTitle(title);
//                h_RecHit_layer[nlayers][eee]->SetMinimum(-10.);
//                h_RecHit_layer[nlayers][eee]->SetMaximum(10.);
			if(eee < 1000) {
				for(int iv = -7; iv < 8; iv++) {
					for(int iu = -7; iu < 8; iu++) {
						if(!IsCellValid.iu_iv_valid(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
//Some thought needs to be put in about the binning and limits of this 1D histogram, probably different for beam type Fermilab and cern.
						CellXY = TheCell.GetCellCoordinatesForPlots(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
						int NumberOfCellVertices = CellXY.size();
						iii = 0;
						if(NumberOfCellVertices == 4) {
							for(it = CellXY.begin(); it != CellXY.end(); it++) {
								HalfHexX[iii] =  it->first;
								HalfHexY[iii++] =  it->second;
							}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
							h_RecHit_layer[nlayers][eee]->AddBin(NumberOfCellVertices, HalfHexX, HalfHexY);
						} else if(NumberOfCellVertices == 6) {
							iii = 0;
							for(it = CellXY.begin(); it != CellXY.end(); it++) {
								FullHexX[iii] =  it->first;
								FullHexY[iii++] =  it->second;
							}
							h_RecHit_layer[nlayers][eee]->AddBin(NumberOfCellVertices, FullHexX, FullHexY);
						}

					}
				}
			}//if eee < 128
		}// loop over events ends here
	}//loop over layers end here


}//contructor ends here


RecHitPlotter_HighGain_New::~RecHitPlotter_HighGain_New()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
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
	edm::Handle<HGCalTBTrackCollection> Tracks;
	event.getByToken(HGCalTBTrackCollection_, Tracks);

	double Average_Pedestal_Per_Event1_Full = 0, Average_Pedestal_Per_Event1_Half = 0, Average_Pedestal_Per_Event2_Full = 0, Average_Pedestal_Per_Event2_Half = 0;
	int Cell_counter1_Full = 0, Cell_counter2_Full = 0, Cell_counter1_Half = 0, Cell_counter2_Half = 0;
	double iux_Max = 0., iyy_Max = 0.;
	int ADC_TMP = 0;

	int Sum_Cluster_Tmp = 0;
//        double iux_CG=0., iyy_CG=0.;
//        int CG_counter=0;
	if(event.id().event() == 1) {
		Time_tmp = event.time().value();
	} else {
		Time_int = event.time().value() - Time_tmp;
		Time_tmp = event.time().value();
//              cout<<endl<<" Event = "<<event.id().event()<<" Time = "<<Time_int<<endl;
		Time_hist->Fill(Time_int);
	}


	for(auto Track : *Tracks) {
		cout << endl << " Event/Trigger = " << event.id().event() << " Track intercept X= " << Track.vertex().X() << " Track intercept Y= " << Track.vertex().Y() << " Slope X= " << Track.momentum().X() << " Slope Y= " << Track.momentum().Y() << " Track hit X= " << Track.pointAt(ExtrapolateZ).X() << " Track hit Y= " << Track.pointAt(ExtrapolateZ).Y() << endl;
	}

	for(auto RecHit1 : *Rechits1) {
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit1.id()).layer(), (RecHit1.id()).sensorIU(), (RecHit1.id()).sensorIV(), (RecHit1.id()).iu(), (RecHit1.id()).iv(), sensorsize);
		double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
		double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);

//             if((RecHit1.energyHigh() > ADC_TMP) && ((RecHit1.id()).cellType() == 0) ){
		if((RecHit1.energyHigh() > ADC_TMP)) {

			ADC_TMP = RecHit1.energyHigh();
			iux_Max = iux;
			iyy_Max = iyy;
		}

//             if(UP && fabs(iux - iux_Max) > 0. && fabs(iyy - iyy_Max) > 0. && RecHit1.energyHigh()< 2000000.){
		if(((RecHit1.id()).cellType() == 0) || ((RecHit1.id()).cellType() == 5)) {

			if((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5)) {
				Cell_counter1_Full++;
				Average_Pedestal_Per_Event1_Full += RecHit1.energyHigh();
			}

			if((iux >= -0.25 && iyy <= -0.5) || (iux >= 0.5)) {
				Cell_counter2_Full++;
				Average_Pedestal_Per_Event2_Full += RecHit1.energyHigh();
			}
		}

		if(((RecHit1.id()).cellType() == 2) || ((RecHit1.id()).cellType() == 3) || ((RecHit1.id()).cellType() == 4) ) {

			if(((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5))) {
				Cell_counter1_Half++;
				Average_Pedestal_Per_Event1_Half += RecHit1.energyHigh();
			}
			if(((iux >= -0.25 && iyy <= -0.5) || (iux >= 0.5))) {
				Cell_counter2_Half++;
				Average_Pedestal_Per_Event2_Half += RecHit1.energyHigh();
			}
		}

//               }


		if(!UP) {
			if(((RecHit1.id()).cellType() == 0)) {

				if((iux <= -0.25 && iux >= -3.25) && (iyy <= -0.25 && iyy >= -5.25 ) && (RecHit1.id()).cellType() == 0) {
					Cell_counter1_Full++;
					Average_Pedestal_Per_Event1_Full += RecHit1.energyHigh();
				}

				if((iux >= -0.25 && iux <= 3.25) && (iyy <= -0.25 && iyy >= -5.25 ) && (RecHit1.id()).cellType() == 0) {
					Cell_counter2_Full++;
					Average_Pedestal_Per_Event2_Full += RecHit1.energyHigh();
				}
			}

			if(((RecHit1.id()).cellType() == 2) || ((RecHit1.id()).cellType() == 3) || ((RecHit1.id()).cellType() == 4) ) {

				if((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5) ) {
					Cell_counter1_Half++;
					Average_Pedestal_Per_Event1_Half += RecHit1.energyHigh();
				}
				if((iux >= 0.25 && iyy <= -0.5) || (iux >= 0.5) ) {
					Cell_counter2_Half++;
					Average_Pedestal_Per_Event2_Half += RecHit1.energyHigh();
				}
			}

		}


	}
	if(ADC_TMP > 5.) {
		CG_X->Fill(iux_Max);
		CG_Y->Fill(iyy_Max);
//   Sum_Cluster_Max->Fill(ADC_TMP);
//   cout<<endl<<" X= "<<CG_X<<" Y= "<<CG_Y<<" Max ADC= "<<ADC_TMP<<endl;
	}

	if(Time_int < 600) ofs << event.id().event() << " " << iux_Max << " " << iyy_Max << endl;


//cout<<endl<<" Event = "<<event.id().event()<<" X max = "<<iux_Max<<" Y max = "<<iyy_Max<<" Energy Max = "<<ADC_TMP<<" FH 1= "<<Cell_counter1_Full<<" FH 2= "<<Cell_counter2_Full<<" HH 1= "<<Cell_counter1_Half<<" HH 2= "<<Cell_counter2_Half<<endl;

	for(auto RecHit : *Rechits) {
		if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
		int n_layer = (RecHit.id()).layer();
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize);
		double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
		double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
		uint32_t EID = essource_.emap_.detId2eid(RecHit.id());
		HGCalTBElectronicsId eid(EID);
		AllCells_Ped->Fill(RecHit.energyHigh());
		if(!DoCommonMode) {
			ADC_Chan[64 * (eid.iskiroc() - 1) + eid.ichan()][event.id().event()] = RecHit.energyHigh();
			h_RecHit_layer[n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, RecHit.energyHigh());
		}
		if(((RecHit.id()).cellType() == 0 ) || ((RecHit.id()).cellType() == 5) ) {
//                  if((RecHit.id()).cellType() == 0){
			if((iux < (iux_Max + 2)  && iux > (iux_Max - 2)) && (iyy < (iyy_Max + 2) && iyy > (iyy_Max - 2)) && RecHit.energyHigh() > 5.  ) {
				Sum_Cluster_Tmp += RecHit.energyHigh();
			}
			if((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5) ) {
				if(!PED && DoCommonMode) {
					if((RecHit.id()).cellType() == 0)  h_RecHit_layer[n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Full / (Cell_counter1_Full))) );
					else if((RecHit.id()).cellType() == 5) h_RecHit_layer[n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Full / (Cell_counter1_Full)) - (Average_Pedestal_Per_Event1_Half / (Cell_counter1_Half))) );
//                                      else if((RecHit.id()).cellType() == 5) h_RecHit_layer[n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (RecHit.energyHigh() - 1.5*(Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full)) ));
				}
				h_RecHit_layer_summed[n_layer - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Full / (Cell_counter1_Full)));
				if(DoCommonMode) {
					ADC_Chan[64 * (eid.iskiroc() - 1) + eid.ichan()][event.id().event()] = RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Full / (Cell_counter1_Full));
					AllCells_CM->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Full / (Cell_counter1_Full)));
				}
//                          Sum_Cluster_ADC->Fill(RecHit.energyHigh()- (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full)));
//                          CG_X->Fill(iux);
//                          CG_Y->Fill(iyy);
			} else if((iux > -0.25 && iyy < -0.50 ) || (iux > 0.50)) {
				if(!PED && DoCommonMode) {
					if((RecHit.id()).cellType() == 0) h_RecHit_layer[n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Full / (Cell_counter2_Full))) );
					else if((RecHit.id()).cellType() == 5) h_RecHit_layer[n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Full / (Cell_counter2_Full)) - (Average_Pedestal_Per_Event2_Half / (Cell_counter2_Half)) ) );

				}
				h_RecHit_layer_summed[n_layer - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Full / (Cell_counter2_Full)));
				if(DoCommonMode) {
					ADC_Chan[64 * (eid.iskiroc() - 1) + eid.ichan()][event.id().event()] = RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Full / (Cell_counter2_Full));
					AllCells_CM->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Full / (Cell_counter2_Full)));
				}
//                          Sum_Cluster_ADC->Fill(RecHit.energyHigh()- (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full)));
//                          CG_X->Fill(iux);
//                          CG_Y->Fill(iyy);
			}
		}
		if((RecHit.id()).cellType() == 1 ) {
			if(!PED) h_RecHit_layer[n_layer - 1][event.id().event()]->Fill(iux , iyy, RecHit.energyHigh());
			if(DoCommonMode) {
				ADC_Chan[64 * (eid.iskiroc() - 1) + eid.ichan()][event.id().event()] = RecHit.energyHigh();
				AllCells_CM->Fill(RecHit.energyHigh());
			}
		}
		if(((RecHit.id()).cellType() != 5) && ((RecHit.id()).cellType() != 1) && ((RecHit.id()).cellType() != 0)) {
			if((iux <= 0.25 && iyy >= -0.25 ) || (iux < -0.5) ) {
				if(!PED && DoCommonMode) h_RecHit_layer[n_layer - 1][event.id().event()]->Fill(iux , iyy, RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Half / Cell_counter1_Half) );
				h_RecHit_layer_summed[n_layer - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Half / Cell_counter1_Half));
				if(DoCommonMode) {
					ADC_Chan[64 * (eid.iskiroc() - 1) + eid.ichan()][event.id().event()] = RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Half / (Cell_counter1_Half));
					AllCells_CM->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event1_Half / (Cell_counter1_Half)));
				}
			}
			if((iux > -0.25 && iyy < -0.50 ) || (iux > 0.50) ) {
				if(!PED && DoCommonMode) h_RecHit_layer[n_layer - 1][event.id().event()]->Fill(iux , iyy, RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Half / Cell_counter2_Half) );
				h_RecHit_layer_summed[n_layer - 1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Half / Cell_counter2_Half));
				if(DoCommonMode) {
					ADC_Chan[64 * (eid.iskiroc() - 1) + eid.ichan()][event.id().event()] = RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Half / (Cell_counter2_Half));
					AllCells_CM->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event2_Half / (Cell_counter2_Half)));
				}
			}

		}

	}

	if(Sum_Cluster_Tmp > 24. && Time_int < 800) Sum_Cluster_ADC->Fill(Sum_Cluster_Tmp);


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
	if(PED) {
		for(int Ievent = 1; Ievent <= 8192; Ievent++) {
			for(int i = 0; i < 128; i++) {
				Mean_SUM[i] += ADC_Chan[i][Ievent - 1];
				Mean_SQ_SUM[i] += pow(ADC_Chan[i][Ievent - 1], 2);
				for(int j = 0; j < 128; j++) {
					Mean_PROD_SUM[i][j] +=  (ADC_Chan[i][Ievent - 1]) * (ADC_Chan[j][Ievent - 1]);
					Diff_IJ_SUM[i][j] += (ADC_Chan[i][Ievent - 1]) - (ADC_Chan[j][Ievent - 1]);
				}
//             cout<<endl<<" ADC Channel i= "<<ADC_Chan[i][Ievent-1]<<endl;
			}
		}

		for(int i = 0; i < 128; i++) {
			Mean_i = Mean_SUM[i] / 8192;
			RMS_i = Return_RMS(Mean_SQ_SUM[i] / 8192, Mean_i);
			for(int j = 0; j < 128; j++) {
				HGCalTBElectronicsId eid2((j / 64) + 1, j - (j / 64) * 64);
				HGCalTBDetId did2 = essource_.emap_.eid2detId(eid2);
				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((did2).layer(), (did2).sensorIU(), (did2).sensorIV(), (did2).iu(), (did2).iv(), sensorsize);
				double iux2 = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
				double iyy2 = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
				Mean_j = Mean_SUM[j] / 8192;
				RMS_j = Return_RMS(Mean_SQ_SUM[j] / 8192, Mean_j);
				Covariance = ((Mean_PROD_SUM[i][j] / 8192) - Mean_i * Mean_j);
				Correlation = ((Mean_PROD_SUM[i][j] / 8192) - Mean_i * Mean_j) / (RMS_i * RMS_j);
				Covar_hist->SetBinContent(i, j, Covariance);
				if(i != j) Correl_hist->SetBinContent(i, j, Correlation);
				h_RecHit_layer[0][i]->Fill(iux2 , iyy2, Correlation);
//            h_RecHit_layer[0][i]->Fill(iux2 ,iyy2,Covariance);
			}
		}
	}

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
