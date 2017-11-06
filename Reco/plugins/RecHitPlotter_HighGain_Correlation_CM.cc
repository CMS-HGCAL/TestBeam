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
  int Eventnumber=0;
  int Eventnumber_1=Eventnumber;
  int pre_Layer=0;
  /*double h_p_Channel[128][6000][MAXLAYERS];
  double l_p_Channel[128][6000][MAXLAYERS];
  double h_A_Channel[128][6000][MAXLAYERS];
  double l_A_Channel[128][6000][MAXLAYERS];*/
  double ADC_Channel[4][128][6000][MAXLAYERS];
  TH2F* h_pre_CM_Final[MAXLAYERS];
  TH2F* l_pre_CM_Final[MAXLAYERS];
  TH2F* h_post_CM_Final[MAXLAYERS];
  TH2F* l_post_CM_Final[MAXLAYERS];
  TH2F* l_pre_CM_Final_Compare[MAXLAYERS];
  TH2F* h_pre_CM_AllLayer_Final;
  TH2F* l_pre_CM_AllLayer_Final;
  TH2F* h_post_CM_AllLayer_Final;
  TH2F* l_post_CM_AllLayer_Final;
  //TH2F* Calculate_Correlation1= new TH2F("Calculate_Correlation1","correlation",2000, -1000., 1000., 2000, -1000., 1000.);
  //TH2F* Calculate_Correlation2= new TH2F("Calculate_Correlation2","correlation",2000, -1000., 1000., 2000, -1000., 1000.);
  //TH2F* Calculate_Correlation3= new TH2F("Calculate_Correlation3","correlation",2000, -1000., 1000., 2000, -1000., 1000.);
  //TH2F* Calculate_Correlation4= new TH2F("Calculate_Correlation4","correlation",2000, -1000., 1000., 2000, -1000., 1000.);
	TH1F* h_Full_Cell[2][MAXLAYERS];
	TH1F* h_Half_Cell[2][MAXLAYERS];
	TH1F* h_MB_Cell[2][MAXLAYERS];
	TH1F* h_Calib_Pads[2][MAXLAYERS];
	TH1F* h_Merged_Cell[2][MAXLAYERS];
  TH1F* l_Full_Cell[2][MAXLAYERS];
  TH1F* l_Half_Cell[2][MAXLAYERS];
  TH1F* l_MB_Cell[2][MAXLAYERS];
  TH1F* l_Calib_Pads[2][MAXLAYERS];
  TH1F* l_Merged_Cell[2][MAXLAYERS];
	TH1F  *h_digi_layer_channel[MAXSKIROCS][64][MAXLAYERS];
//        TH2F  *h_digi_layer_channel_CM[2][64];
	TH1F* Sum_Cluster_ADC;
	TH1F* AllCells_Ped;
	TH1F* AllCells_CM;
  TH1F* h_Correlation_CellArea[2][MAXLAYERS];
  TH1F* h_Correlation_CellPerimeter[2][MAXLAYERS];
  TH1F* l_Correlation_CellArea[2][MAXLAYERS];
  TH1F* l_Correlation_CellPerimeter[2][MAXLAYERS];
	TH2F* Noise_2D_Profile;
  char name[50], title[50],xtitle[50],ytitle[50];
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
	  for(int ISkiroc=1;ISkiroc<3;++ISkiroc){
	    //printf("ISkiroc=%d\n",ISkiroc);
	    sprintf(name, "HighGainCommonModeVsCellArea_%i", ILayer+1);
	    sprintf(title, "HighGainCommonMode Noise Vs CellArea Layer %i", ILayer+1);
	    h_Correlation_CellArea[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, 0., 1.5);
	    h_Correlation_CellArea[ISkiroc-1][ILayer]->GetXaxis()->SetTitle("Area[cm^{2}]");
	    h_Correlation_CellArea[ISkiroc-1][ILayer]->GetYaxis()->SetTitle("Common Mode Noise");
	    h_Correlation_CellArea[ISkiroc-1][ILayer]->SetMarkerStyle(20);

	    sprintf(name, "HighGainCommonModeVsCellPerimeter_%i", ILayer+1);
	    sprintf(title, "HighGainCommonMode Noise Vs CellPerimeter Layer %i", ILayer+1);
	    h_Correlation_CellPerimeter[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, 0., 5.);
	    h_Correlation_CellPerimeter[ISkiroc-1][ILayer]->GetXaxis()->SetTitle("Perimeter[cm]");
	    h_Correlation_CellPerimeter[ISkiroc-1][ILayer]->GetYaxis()->SetTitle("Common Mode Noise");
	    h_Correlation_CellPerimeter[ISkiroc-1][ILayer]->SetMarkerStyle(20);
	    
	    sprintf(name, "LowGainCommonModeVsCellArea_%i", ILayer+1);
	    sprintf(title, "LowgainCommonMode Noise Vs CellArea Layer %i", ILayer+1);
	    l_Correlation_CellArea[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, 0., 1.5);
	    l_Correlation_CellArea[ISkiroc-1][ILayer]->GetXaxis()->SetTitle("Area[cm^{2}]");
	    l_Correlation_CellArea[ISkiroc-1][ILayer]->GetYaxis()->SetTitle("Common Mode Noise");
	    l_Correlation_CellArea[ISkiroc-1][ILayer]->SetMarkerStyle(20);

	    sprintf(name, "LowGainCommonModeVsCellPerimeter_%i", ILayer+1);
	    sprintf(title, "LoGainCommonMode Noise Vs CellPerimeter Layer %i", ILayer+1);
	    l_Correlation_CellPerimeter[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, 0., 5.);
	    l_Correlation_CellPerimeter[ISkiroc-1][ILayer]->GetXaxis()->SetTitle("Perimeter[cm]");
	    l_Correlation_CellPerimeter[ISkiroc-1][ILayer]->GetYaxis()->SetTitle("Common Mode Noise");
	    l_Correlation_CellPerimeter[ISkiroc-1][ILayer]->SetMarkerStyle(20);
	  
	    sprintf(name, "H_Full_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
	    sprintf(title, "H_Full Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		h_Full_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "H_Half_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "H_Half Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		h_Half_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "H_MB_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "H_MB Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		h_MB_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "H_Calib_Pads_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "H_Calib Pads Layer %i Ski %d", ILayer+1,ISkiroc);
		h_Calib_Pads[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "H_Merged_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "H_Merged Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		h_Merged_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);

		sprintf(name, "L_Full_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "L_Full Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		l_Full_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "L_Half_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "L_Half Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		l_Half_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "L_MB_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "L_MB Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		l_MB_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "L_Calib_Pads_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "L_Calib Pads Layer %i Ski %d", ILayer+1,ISkiroc);
		l_Calib_Pads[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
		sprintf(name, "L_Merged_Cell_Layer_%i_Ski%d", ILayer+1,ISkiroc);
		sprintf(title, "L_Merged Cell Layer %i Ski %d", ILayer+1,ISkiroc);
		l_Merged_Cell[ISkiroc-1][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
	  }
		for(int ISkiroc = 1; ISkiroc <= MAXSKIROCS; ISkiroc++) {
			for(int Channel = 0; Channel < 64; Channel++) {
				sprintf(name, "Ski_%i_Channel_%i_Layer_%i", ISkiroc, Channel, ILayer+1);
				sprintf(title, "Ski %i Channel %i Layer %i", ISkiroc, Channel, ILayer+1);
				h_digi_layer_channel[ISkiroc - 1][Channel][ILayer] = fs->make<TH1F>(name, title, 1000, -500., 500.);
				/*
				        			 sprintf(name, "Ski_%i_Channel_%i_CM",ISkiroc,Channel);
				       			         sprintf(title, "Ski %i Channel %i CM",ISkiroc,Channel);
				                		 h_digi_layer_channel_CM[ISkiroc-1][Channel] = fs->make<TH2F>(name, title, 1000,-500., 500.,1000,-500., 500.);
				*/
			}
		}
		sprintf(name,"h_pre_CM_correlation_Layer%d",ILayer+1);
		sprintf(title,"h_pre_CM_correlation_Layer%d",ILayer+1);
		h_pre_CM_Final[ILayer] = fs->make<TH2F>(name,title, 128, 0, 128, 128, 0 ,128);
		h_pre_CM_Final[ILayer] -> GetXaxis()-> SetTitle("Channel");
		h_pre_CM_Final[ILayer] -> GetYaxis()-> SetTitle("Channel");

		sprintf(name,"h_post_CM_correlation_Layer%d",ILayer+1);
                sprintf(title,"h_post_CM_correlation_Layer%d",ILayer+1);
                h_post_CM_Final[ILayer] = fs->make<TH2F>(name,title, 128, 0, 128, 128, 0 ,128);
                h_post_CM_Final[ILayer] -> GetXaxis()-> SetTitle("Channel");
                h_post_CM_Final[ILayer] -> GetYaxis()-> SetTitle("Channel");

		sprintf(name,"l_pre_CM_correlation_Layer%d",ILayer+1);
		sprintf(title,"l_pre_CM_correlation_Layer%d",ILayer+1);
                l_pre_CM_Final[ILayer] = fs->make<TH2F>(name,title, 128, 0, 128, 128, 0 ,128);
                l_pre_CM_Final[ILayer] -> GetXaxis()-> SetTitle("Channel");
		l_pre_CM_Final[ILayer] -> GetYaxis()-> SetTitle("Channel");

                sprintf(name,"l_post_CM_correlation_Layer%d",ILayer+1);
                sprintf(title,"l_post_CM_correlation_Layer%d",ILayer+1);
                l_post_CM_Final[ILayer] = fs->make<TH2F>(name,title, 128, 0, 128, 128, 0 ,128);
		l_post_CM_Final[ILayer] -> GetXaxis()-> SetTitle("Channel");
                l_post_CM_Final[ILayer] -> GetYaxis()-> SetTitle("Channel");

		//sprintf(name,"h_post_CM_correlation_Layer%d_C",ILayer+1);
                //sprintf(title,"h_post_CM_correlation_Layer%d_C",ILayer+1);
                //l_pre_CM_Final_Compare[ILayer] = fs->make<TH2F>(name,title, 128, 0, 127, 128, 0 ,127);
                //l_pre_CM_Final_Compare[ILayer] -> GetXaxis()-> SetTitle("Channel");
		//l_pre_CM_Final_Compare[ILayer] -> GetYaxis()-> SetTitle("Channel");
	}
	sprintf(name,"h_pre_CM_correlation_AllLayers");
	sprintf(title,"h_pre_CM_correlation_AllLayers");
	h_pre_CM_AllLayer_Final = fs->make<TH2F>(name,title, 128*8, 0, 128*8, 128*8, 0 ,128*8);
	h_pre_CM_AllLayer_Final -> GetXaxis()-> SetTitle("Channel");
	h_pre_CM_AllLayer_Final -> GetYaxis()-> SetTitle("Channel");

	sprintf(name,"h_post_CM_correlation_AllLayers");
	sprintf(title,"h_post_CM_correlation_AllLayers");
	h_post_CM_AllLayer_Final = fs->make<TH2F>(name,title, 128*8, 0, 128*8, 128*8, 0 ,128*8);
	h_post_CM_AllLayer_Final -> GetXaxis()-> SetTitle("Channel");
	h_post_CM_AllLayer_Final -> GetYaxis()-> SetTitle("Channel");

	sprintf(name,"l_pre_CM_correlation_AllLayers");
	sprintf(title,"l_pre_CM_correlation_AllLayers");
	l_pre_CM_AllLayer_Final = fs->make<TH2F>(name,title, 128*8, 0, 128*8, 128*8, 0 ,128*8);
	l_pre_CM_AllLayer_Final -> GetXaxis()-> SetTitle("Channel");
	l_pre_CM_AllLayer_Final -> GetYaxis()-> SetTitle("Channel");

	sprintf(name,"l_post_CM_correlation_AllLayers");
	sprintf(title,"l_post_CM_correlation_AllLayers");
	l_post_CM_AllLayer_Final = fs->make<TH2F>(name,title, 128*8, 0, 128*8, 128*8, 0 ,128*8);
	l_post_CM_AllLayer_Final -> GetXaxis()-> SetTitle("Channel");
	l_post_CM_AllLayer_Final -> GetYaxis()-> SetTitle("Channel");
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

	double h_Average_Pedestal_Per_Event_Full[2][MAXLAYERS]           = {{0}};
	double l_Average_Pedestal_Per_Event_Full[2][MAXLAYERS]           = {{0}};
	int    Cell_counter[2][MAXLAYERS]                                = {{0}};   
	double h_Average_Pedestal_Per_Event_Half[2][MAXLAYERS]           = {{0}};
	double l_Average_Pedestal_Per_Event_Half[2][MAXLAYERS]           = {{0}};
	int    Cell_counter_Half[2][MAXLAYERS]                           = {{0}};
	double h_Average_Pedestal_Per_Event_MB[2][MAXLAYERS]             = {{0}};
	double l_Average_Pedestal_Per_Event_MB[2][MAXLAYERS]             = {{0}};
	int    Cell_counter_MB[2][MAXLAYERS]                             = {{0}};
	double h_Average_Pedestal_Per_Event_Calib_Pad[2][MAXLAYERS]      = {{0}};
	double l_Average_Pedestal_Per_Event_Calib_Pad[2][MAXLAYERS]      = {{0}};
	int    Cell_counter_Calib_Pad[2][MAXLAYERS]                      = {{0}};
	double h_Average_Pedestal_Per_Event_Merged_Cell[2][MAXLAYERS]    = {{0}};
	double l_Average_Pedestal_Per_Event_Merged_Cell[2][MAXLAYERS]    = {{0}};
	int    Cell_counter_Merged_Cell[2][MAXLAYERS]                    = {{0}};
	for(auto RecHit1 : *Rechits1) {
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit1.id()).layer(), (RecHit1.id()).sensorIU(), (RecHit1.id()).sensorIV(), (RecHit1.id()).iu(), (RecHit1.id()).iv(), sensorsize);
		uint32_t EID = essource_.emap_.detId2eid(RecHit1.id());
		HGCalTBElectronicsId eid(EID);
		if(RecHit1.energyHigh() > 30) continue;
		if((RecHit1.id()).cellType() == 0) {
		  Cell_counter[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1]++;
		  h_Average_Pedestal_Per_Event_Full[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		  l_Average_Pedestal_Per_Event_Full[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyLow();
		}
		else if((RecHit1.id()).cellType() == 2) {
		  Cell_counter_Half[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1]++;
		  h_Average_Pedestal_Per_Event_Half[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		  l_Average_Pedestal_Per_Event_Half[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyLow();
		}
		else if(((RecHit1.id()).cellType() == 3)){
		  Cell_counter_MB[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1]++;
		  h_Average_Pedestal_Per_Event_MB[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		  l_Average_Pedestal_Per_Event_MB[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyLow();
		}
		else if(((RecHit1.id()).cellType() == 5)) {
		  Cell_counter_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1]++;
		  h_Average_Pedestal_Per_Event_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		  l_Average_Pedestal_Per_Event_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyLow();
		}
		else if((RecHit1.id()).cellType() == 1) {
		  Cell_counter_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1]++;
		  h_Average_Pedestal_Per_Event_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyHigh();
		  l_Average_Pedestal_Per_Event_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit1.id()).layer() - 1] += RecHit1.energyLow();
		}

	}


	for(int iii = 0; iii < MAXLAYERS; iii++) {
	  for(int ISkiroc=0;ISkiroc<2;++ISkiroc){
	  if(Cell_counter[ISkiroc][iii] != 0) h_Full_Cell[ISkiroc][iii]->Fill(h_Average_Pedestal_Per_Event_Full[ISkiroc][iii] / Cell_counter[ISkiroc][iii]);
	  if(Cell_counter_Half[ISkiroc][iii] != 0) h_Half_Cell[ISkiroc][iii]->Fill(h_Average_Pedestal_Per_Event_Half[ISkiroc][iii] / Cell_counter_Half[ISkiroc][iii]);
	  if(Cell_counter_MB[ISkiroc][iii] != 0) h_MB_Cell[ISkiroc][iii]->Fill(h_Average_Pedestal_Per_Event_MB[ISkiroc][iii] / Cell_counter_MB[ISkiroc][iii]);
	  if(Cell_counter_Merged_Cell[ISkiroc][iii] != 0) h_Merged_Cell[ISkiroc][iii]->Fill(h_Average_Pedestal_Per_Event_Merged_Cell[ISkiroc][iii] / Cell_counter_Merged_Cell[ISkiroc][iii]);
	  if(Cell_counter_Calib_Pad[ISkiroc][iii] != 0) h_Calib_Pads[ISkiroc][iii]->Fill(h_Average_Pedestal_Per_Event_Calib_Pad[ISkiroc][iii] / Cell_counter_Calib_Pad[ISkiroc][iii]);
	  
	  if(Cell_counter[ISkiroc][iii] != 0) l_Full_Cell[ISkiroc][iii]->Fill(l_Average_Pedestal_Per_Event_Full[ISkiroc][iii] / Cell_counter[ISkiroc][iii]);
	  if(Cell_counter_Half[ISkiroc][iii] != 0) l_Half_Cell[ISkiroc][iii]->Fill(l_Average_Pedestal_Per_Event_Half[ISkiroc][iii] / Cell_counter_Half[ISkiroc][iii]);
	  if(Cell_counter_MB[ISkiroc][iii] != 0) l_MB_Cell[ISkiroc][iii]->Fill(l_Average_Pedestal_Per_Event_MB[ISkiroc][iii] / Cell_counter_MB[ISkiroc][iii]);
	  if(Cell_counter_Merged_Cell[ISkiroc][iii] != 0) l_Merged_Cell[ISkiroc][iii]->Fill(l_Average_Pedestal_Per_Event_Merged_Cell[ISkiroc][iii] / Cell_counter_Merged_Cell[ISkiroc][iii]);
	  if(Cell_counter_Calib_Pad[ISkiroc][iii] != 0) l_Calib_Pads[ISkiroc][iii]->Fill(l_Average_Pedestal_Per_Event_Calib_Pad[ISkiroc][iii] / Cell_counter_Calib_Pad[ISkiroc][iii]);
	  }
	}

	for(auto RecHit : *Rechits) {
		if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize);
		uint32_t EID = essource_.emap_.detId2eid(RecHit.id());
		HGCalTBElectronicsId eid(EID);

		if((RecHit.id()).layer() - 1 != pre_Layer){ 
		  if((RecHit.id()).layer() - 1 != 0) Eventnumber=Eventnumber_1;
		  else Eventnumber_1=Eventnumber;
		}
		
		if(!doCommonMode_CM) {
			h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh());
			Noise_2D_Profile->Fill((64 * (eid.iskiroc() - 1) + eid.ichan()), RecHit.energyHigh());
		}
		if(doCommonMode_CM) {
			AllCells_Ped->Fill(RecHit.energyHigh());
			ADC_Channel[0][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyHigh();
			ADC_Channel[1][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyLow();
			if((RecHit.id()).cellType() == 0 || (RecHit.id()).cellType() == 4){
			  h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Full[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1])));
			  ADC_Channel[2][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Full[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			  ADC_Channel[3][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyLow() - (l_Average_Pedestal_Per_Event_Full[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			}
			else if((RecHit.id()).cellType() == 2 ){
			  h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Half[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Half[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1])));
			  ADC_Channel[2][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Half[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Half[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			  ADC_Channel[3][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyLow() - (l_Average_Pedestal_Per_Event_Half[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Half[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			}
			else if((RecHit.id()).cellType() == 1){
			  h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1])));
			  ADC_Channel[2][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
                          ADC_Channel[3][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyLow() - (l_Average_Pedestal_Per_Event_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Calib_Pad[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			}
			else if(((RecHit.id()).cellType() == 3) ) {
			  h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_MB[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_MB[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1])));
			  ADC_Channel[2][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_MB[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_MB[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			  ADC_Channel[3][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyLow() - (l_Average_Pedestal_Per_Event_MB[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_MB[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			}
			 else if(((RecHit.id()).cellType() == 5)) {
			  h_digi_layer_channel[eid.iskiroc() - 1][eid.ichan()][(RecHit.id()).layer() - 1]->Fill(RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1])));
			  ADC_Channel[2][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			  ADC_Channel[3][(64 *((eid.iskiroc() - 1)%2))+eid.ichan()][Eventnumber][(RecHit.id()).layer() - 1]= RecHit.energyLow() - (l_Average_Pedestal_Per_Event_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter_Merged_Cell[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1]));
			}

			Noise_2D_Profile->Fill((64 * (eid.iskiroc() - 1) + eid.ichan()), RecHit.energyHigh() - (h_Average_Pedestal_Per_Event_Full[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1] / (Cell_counter[(eid.iskiroc() - 1)%2][(RecHit.id()).layer() - 1])));
		}

		pre_Layer=(RecHit.id()).layer() - 1;

	}
	++Eventnumber;
	//printf("analyze method ends\n");
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
  for(int iii = 0; iii < MAXLAYERS; iii++) {
    for(int ISkiroc=0;ISkiroc<2;++ISkiroc){
    h_Correlation_CellArea[ISkiroc][iii]->Fill(1.579, h_Merged_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellArea[ISkiroc][iii]->Fill(1.09,  h_Full_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellArea[ISkiroc][iii]->Fill(0.978, h_MB_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellArea[ISkiroc][iii]->Fill(0.545, h_Half_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellArea[ISkiroc][iii]->Fill(0.15,  h_Calib_Pads[ISkiroc][iii]->GetRMS());

    h_Correlation_CellPerimeter[ISkiroc][iii]->Fill(5.972, h_Merged_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellPerimeter[ISkiroc][iii]->Fill(3.894, h_Full_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellPerimeter[ISkiroc][iii]->Fill(3.47,  h_Half_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellPerimeter[ISkiroc][iii]->Fill(4.156, h_MB_Cell[ISkiroc][iii]->GetRMS());
    h_Correlation_CellPerimeter[ISkiroc][iii]->Fill(1.44,  h_Calib_Pads[ISkiroc][iii]->GetRMS());
    
    l_Correlation_CellArea[ISkiroc][iii]->Fill(1.579, l_Merged_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellArea[ISkiroc][iii]->Fill(1.09,  l_Full_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellArea[ISkiroc][iii]->Fill(0.978, l_MB_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellArea[ISkiroc][iii]->Fill(0.545, l_Half_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellArea[ISkiroc][iii]->Fill(0.15,  l_Calib_Pads[ISkiroc][iii]->GetRMS());

    l_Correlation_CellPerimeter[ISkiroc][iii]->Fill(5.972, l_Merged_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellPerimeter[ISkiroc][iii]->Fill(3.894, l_Full_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellPerimeter[ISkiroc][iii]->Fill(3.47,  l_Half_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellPerimeter[ISkiroc][iii]->Fill(4.156, l_MB_Cell[ISkiroc][iii]->GetRMS());
    l_Correlation_CellPerimeter[ISkiroc][iii]->Fill(1.44,  l_Calib_Pads[ISkiroc][iii]->GetRMS());
    }
  }

  printf("Eventnumber=%d\n",Eventnumber);

  /* //=========================One Layer Correlation========================
  for(int layer=0;layer<8;layer++){
    if(Eventnumber==0)break; 
    for(int ChX=0;ChX<128;ChX++){
      // Calculate_Correlation-> Reset();
      for(int ChY=ChX;ChY<128;ChY++){
	double sumxy=0,sumx=0,sumy=0,sumx2=0,sumy2=0;
	double correlation=0;
	for(int Event=0;Event<Eventnumber;Event++){
	  //Calculate_Correlation1->Fill(h_p_Channel[ChX][Event][layer],h_p_Channel[ChY][Event][layer]);
	  //Calculate_Correlation2->Fill(l_p_Channel[ChX][Event][layer],l_p_Channel[ChY][Event][layer]);
	  Calculate_Correlation3->Fill(h_A_Channel[ChX][Event][layer],h_A_Channel[ChY][Event][layer]);
          //Calculate_Correlation4->Fill(l_A_Channel[ChX][Event][layer],l_A_Channel[ChY][Event][layer]);

	  sumxy += h_A_Channel[ChX][Event][layer]*h_A_Channel[ChY][Event][layer];
	  sumx += h_A_Channel[ChX][Event][layer];
	  sumy += h_A_Channel[ChY][Event][layer];
	  sumx2 += pow(h_A_Channel[ChX][Event][layer],2);
	  sumy2 += pow(h_A_Channel[ChY][Event][layer],2);
	}
	correlation=((Eventnumber*sumxy)-(sumx*sumy))/(sqrt((Eventnumber*sumx2)-(sumx*sumx))*sqrt((Eventnumber*sumy2)-(sumy*sumy)));
	
	//h_pre_CM_Final[layer]->Fill(ChX,ChY,Calculate_Correlation1-> GetCorrelationFactor());
	//l_pre_CM_Final[layer]->Fill(ChX,ChY,Calculate_Correlation2-> GetCorrelationFactor());
	l_pre_CM_Final_Compare[layer]->Fill(ChX,ChY,correlation);
	h_post_CM_Final[layer]->Fill(ChX,ChY,Calculate_Correlation3-> GetCorrelationFactor());
        //l_post_CM_Final[layer]->Fill(ChX,ChY,Calculate_Correlation4-> GetCorrelationFactor());
	if(ChX!=ChY){
	  //h_pre_CM_Final[layer]->Fill(ChY,ChX,Calculate_Correlation1-> GetCorrelationFactor());
	  //l_pre_CM_Final[layer]->Fill(ChY,ChX,Calculate_Correlation2-> GetCorrelationFactor());
	  l_pre_CM_Final_Compare[layer]->Fill(ChY,ChX,correlation);
	  h_post_CM_Final[layer]->Fill(ChY,ChX,Calculate_Correlation3-> GetCorrelationFactor());
	  //l_post_CM_Final[layer]->Fill(ChY,ChX,Calculate_Correlation4-> GetCorrelationFactor());
	}
	//printf("[%d][%d][%d]=%lf\n",ChX,ChY,layer,Calculate_Correlation1-> GetCorrelationFactor());
	//Calculate_Correlation1-> Reset();
	//Calculate_Correlation2-> Reset();
	Calculate_Correlation3-> Reset();
        //Calculate_Correlation4-> Reset();
      }
      printf("Finish Channel[%d] Layer[%d]\n",ChX,layer+1);
    }
  }
  //==========================================================================*/

  //=========================All Layers Correlation===========================
  for(int layer=0;layer<8;++layer){
    if(Eventnumber==0)break; 
    for(int ChX=0;ChX<128;++ChX){
      // Calculate_Correlation-> Reset();
      for(int layer2=layer;layer2<8;++layer2){
	for(int ChY=ChX;ChY<128;ChY++){
	  double sumxy[4]={0},sumx[4]={0},sumy[4]={0},sumx2[4]={0},sumy2[4]={0},correlation[4]={0};
	  for(int Event=0;Event<Eventnumber;++Event){
	    //Calculate_Correlation1->Fill(h_p_Channel[ChX][Event][layer],h_p_Channel[ChY][Event][layer2]);
	    //Calculate_Correlation2->Fill(l_p_Channel[ChX][Event][layer],l_p_Channel[ChY][Event][layer2]);
	    //Calculate_Correlation3->Fill(h_A_Channel[ChX][Event][layer],h_A_Channel[ChY][Event][layer2]);
	    //Calculate_Correlation4->Fill(l_A_Channel[ChX][Event][layer],l_A_Channel[ChY][Event][layer2]);
	    for(int type=0;type<4;++type){
	      sumxy[type] += ADC_Channel[type][ChX][Event][layer]*ADC_Channel[type][ChY][Event][layer2];
	      sumx[type] += ADC_Channel[type][ChX][Event][layer];
	      sumy[type] += ADC_Channel[type][ChY][Event][layer2];
	      sumx2[type] += pow(ADC_Channel[type][ChX][Event][layer],2);
	      sumy2[type] += pow(ADC_Channel[type][ChY][Event][layer2],2);
	    }
	  }
	  for(int type=0;type<4;++type){
	    correlation[type]=((Eventnumber*sumxy[type])-(sumx[type]*sumy[type]))/(sqrt((Eventnumber*sumx2[type])-(sumx[type]*sumx[type]))*sqrt((Eventnumber*sumy2[type])-(sumy[type]*sumy[type])));
	  }
	  if(layer2==layer){
	    h_pre_CM_Final[layer]->Fill(ChX,ChY,correlation[0]);
	    l_pre_CM_Final[layer]->Fill(ChX,ChY,correlation[1]);
	    h_post_CM_Final[layer]->Fill(ChX,ChY,correlation[2]);
	    l_post_CM_Final[layer]->Fill(ChX,ChY,correlation[3]);
	    if(ChX!=ChY){
	      h_pre_CM_Final[layer]->Fill(ChY,ChX,correlation[0]);
	      l_pre_CM_Final[layer]->Fill(ChY,ChX,correlation[1]);
	      h_post_CM_Final[layer]->Fill(ChY,ChX,correlation[2]);
	      l_post_CM_Final[layer]->Fill(ChY,ChX,correlation[3]);
	    }
	  }
	  h_pre_CM_AllLayer_Final->Fill(ChX+(layer*128),ChY+(layer2*128),correlation[0]);
	  l_pre_CM_AllLayer_Final->Fill(ChX+(layer*128),ChY+(layer2*128),correlation[1]);
	  h_post_CM_AllLayer_Final->Fill(ChX+(layer*128),ChY+(layer2*128),correlation[2]);
	  l_post_CM_AllLayer_Final->Fill(ChX+(layer*128),ChY+(layer2*128),correlation[3]);
	  if(ChX!=ChY){
	    h_pre_CM_AllLayer_Final->Fill(ChY+(layer*128),ChX+(layer2*128),correlation[0]);
	    l_pre_CM_AllLayer_Final->Fill(ChY+(layer*128),ChX+(layer2*128),correlation[1]);
	    h_post_CM_AllLayer_Final->Fill(ChY+(layer*128),ChX+(layer2*128),correlation[2]);
	    l_post_CM_AllLayer_Final->Fill(ChY+(layer*128),ChX+(layer2*128),correlation[3]);
	  }
	  if(layer2!=layer){
	    h_pre_CM_AllLayer_Final->Fill(ChY+(layer2*128),ChX+(layer*128),correlation[0]);
	    l_pre_CM_AllLayer_Final->Fill(ChY+(layer2*128),ChX+(layer*128),correlation[1]);
	    h_post_CM_AllLayer_Final->Fill(ChY+(layer2*128),ChX+(layer*128),correlation[2]);
	    l_post_CM_AllLayer_Final->Fill(ChY+(layer2*128),ChX+(layer*128),correlation[3]);
	    if(ChX!=ChY){
	      h_pre_CM_AllLayer_Final->Fill(ChX+(layer2*128),ChY+(layer*128),correlation[0]);
	      l_pre_CM_AllLayer_Final->Fill(ChX+(layer2*128),ChY+(layer*128),correlation[1]);
	      h_post_CM_AllLayer_Final->Fill(ChX+(layer2*128),ChY+(layer*128),correlation[2]);
	      l_post_CM_AllLayer_Final->Fill(ChX+(layer2*128),ChY+(layer*128),correlation[3]);
	    }
	  }
	}
      }
      printf("Finish Channel[%d] Layer[%d]\n",ChX,layer+1);
    }
  }
  //==================================================================================================
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
