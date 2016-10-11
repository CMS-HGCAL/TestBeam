/* Need full layer, 7 cell cluster, and 19 cell cluster histograms for each layer
 * also need each for all layers summed
 * use ADC to MIP conversion of 1 MIP = 10 ADC Counts */

/**
	@Author: Ryan Quinn <ryan>
		7 July 2016
		quinn@physics.umn.edu
*/



// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include <sstream>
#include <fstream>
#include <math.h>
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
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

// chooses which particle to look at. Inverts threshold filtering.
// if nothing is selected, electrons are the default
const bool ELECTRONS(1);// uses > *CELLS_THRESHOLD
const bool PIONS(0);// uses > PION_*CELLS_THRESHOLD and < *CELLS_THRESHOLD

double Layer_Z[16]  = {1.2, 2., 3.5, 4.3, 5.8, 6.3, 8.7, 9.5, 11.4, 12.2, 13.8, 14.6, 16.6, 17.4, 20., 20.8};

double ADCtoMIP[16] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. };
double MIP2ParticleCalib = 1.3;  // FIXME for CERN

//double ADCtoMIP[16] = {16.02,16.85,15.36,14.73,10.66,15.64,16.52,14.24,10.07,14.42,16.14,17.33,16.61,16.84,15.79,16.43};// one MIP is equal to _ADCtoMIP_ ADC Counts
//double LayerWeight[16] = {0.6374029601923652, 0.7392021202456731, 0.6374273268336504, 0.7392021202456731, 0.6374273268336504, 0.8861075434658853, 0.8487578715427883, 1.0330129666860974, 0.8487578715427883, 1.0330129666860974, 0.8487578715427883, 1.5226977107534714, 1.2714189609610644, 1.5226977107534714, 1.2714189609610644, 1.5226977107534714};// X0 weights

//double LayerWeight[16] = {1.4091566745180932, 0.7020676448403224, 0.6054055986179145, 0.7020676448403224, 0.6054055986179145, 0.8415931435769973, 0.8061197656138868, 0.9811186423136724, 0.8061197656138868, 0.9811186423136724, 0.8061197656138868, 1.4462036381025891, 1.2075480996058319, 1.4462036381025891, 1.2075480996058319, 1.4462036381025891};


double LayerWeight_8L_conf1[16] = {33.074, 13.184, 14.17, 9.788, 9.766, 9.766, 16.339, 14.129};
double X0depth_8L_conf1[16] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42};
double LayerWeight_8L_conf2[16] = {35.866, 30.864, 28.803, 23.095, 20.657, 19.804, 36.322, 27.451};
double X0depth_8L_conf2[16] = {5.048, 3.412, 3.412, 2.866, 2.512, 1.625, 2.368, 6.021};
double weights2GeV = 1.e-03;

//double LayerWeight[16] = {0.4847555727337982, 1.0214605968539232, 0.4847555727337982, 1.0214605968539232, 0.4847555727337982, 1.1420105918768606, 0.6423912113800805, 1.2625605868997982, 0.6423912113800805, 1.2625605868997982, 0.6423912113800805, 1.6643939036429232, 0.9576624886726451, 1.6643939036429232, 0.9576624886726451, 1.6643939036429232};// dE/dx weights

double LayerSumWeight = 1.;
const int CMTHRESHOLD = 30;// anything less than this value is added to the commonmode sum

// applied to all layers sum after commonmode subtraction and the ADC to MIP conversion
const double ALLCELLS_THRESHOLD = 50.;
const double SEVENCELLS_THRESHOLD = 50.;
const double NINETEENCELLS_THRESHOLD = 50.;
const double PION_ALLCELLS_THRESHOLD = 15.;
const double PION_7CELLS_THRESHOLD = -100.;
const double PION_19CELLS_THRESHOLD = -100.;


using namespace std;



class Layer_Sum_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit Layer_Sum_Analyzer(const edm::ParameterSet&);
	~Layer_Sum_Analyzer();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;

	// ----------member data ---------------------------
	edm::EDGetToken HGCalTBRecHitCollection_;
        int CERN_8layers_config_;

	struct {
		HGCalElectronicsMap emap_;
	} essource_;
	string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
	int sensorsize = 128;
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
        HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
        double maxdist = (1 + sqrt (3) / 2) * HGCAL_TB_CELL::FULL_CELL_SIDE;  // <<< FIXME maxdist > HGCAL_TB_CELL::FULL_CELL_SIDE !!

        double Weights_8L[MAXLAYERS];
        double X0_8L[MAXLAYERS];

        TH1F *h_CM_layer[MAXSKIROCS];
	TH1F *h_sum_layer[MAXLAYERS], *h_layer_seven[MAXLAYERS], *h_layer_nineteen[MAXLAYERS], *h_Seed_layer[MAXLAYERS];
	TH1F *h_sum_layer_AbsW[MAXLAYERS], *h_layer_seven_AbsW[MAXLAYERS], *h_layer_nineteen_AbsW[MAXLAYERS], *h_Seed_layer_AbsW[MAXLAYERS];
      	TH1F *h_x_layer[MAXLAYERS], *h_y_layer[MAXLAYERS];
	TH2F *h_x_y_layer[MAXLAYERS];

        TH1F* h_Radius[MAXLAYERS];      
        TH1F *h_E1oE7_layer[MAXLAYERS], *h_E1oE19_layer[MAXLAYERS], *h_E7oE19_layer[MAXLAYERS];

        TH1F *h_sum_all, *h_seven_all, *h_nineteen_all;
        TH1F *h_sum_all_AbsW, *h_seven_all_AbsW, *h_nineteen_all_AbsW;
        TProfile *tp_E1_vs_layer, *tp_E7_vs_layer, *tp_E19_vs_layer;
        TProfile *tp_E1oSumL_vs_layer, *tp_E7oSumL_vs_layer, *tp_E19oSumL_vs_layer;
        TProfile *tp_E1oE7_vs_layer, *tp_E1oE19_vs_layer, *tp_E7oE19_vs_layer;

	TH2F *HighGain_LowGain_2D;
	TH2F *Energy_LowGain_2D;

        int EVENT = 0;

	map<int, double> Time_Stamp;
	map<int, double> Delta_Time_Stamp;
	double Time_Temp = 0.;
};



Layer_Sum_Analyzer::Layer_Sum_Analyzer(const edm::ParameterSet& iConfig)
{
  std::cout << " welcome costruttore MAXSKIROCS = " << MAXSKIROCS << std::endl;
	// initialization
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	CERN_8layers_config_ = iConfig.getParameter<int>("CERN_8layers_config");

	//booking the histos
	for(int layer = 0; layer < MAXLAYERS; layer++) {
	  h_sum_layer[layer] = fs->make<TH1F>(Form("sumAll_Layer%d", layer+1), "", 40010, -10, 40000);
	  h_Seed_layer[layer] = fs->make<TH1F>(Form("h_Seed_layer%d",layer+1), Form("h_Seed_layer%d",layer+1), 40010, -10, 40000);
	  h_layer_seven[layer] = fs->make<TH1F>(Form("sum7_Layer%d", layer+1),"",  40010, -10, 40000);
	  h_layer_nineteen[layer] = fs->make<TH1F>(Form("sum19_Layer%d", layer+1), "", 40010, -10, 40000);

	  h_sum_layer_AbsW[layer] = fs->make<TH1F>(Form("sumAll_Layer%d_AbsW", layer+1), "", 40010, -10, 40000);
	  h_Seed_layer_AbsW[layer] = fs->make<TH1F>(Form("h_Seed_layer%d_AbsW",layer+1), Form("h_Seed_layer%d",layer+1), 40010, -10, 40000);
	  h_layer_seven_AbsW[layer] = fs->make<TH1F>(Form("sum7_Layer%d_AbsW", layer+1),"",  40010, -10, 40000);
	  h_layer_nineteen_AbsW[layer] = fs->make<TH1F>(Form("sum19_Layer%d_AbsW", layer+1), "", 40010, -10, 40000);

	  h_x_layer[layer] = fs->make<TH1F>(Form("X_Layer%d", layer+1), "", 2000, -10., 10. );
	  h_y_layer[layer] = fs->make<TH1F>(Form("Y_Layer%d", layer+1), "", 2000, -10., 10. );
	  h_x_y_layer[layer] = fs->make<TH2F>(Form("YvsX_Layer%d", layer+1), "", 2000, -10., 10., 2000, -10., 10. );

	  h_E1oE7_layer[layer] = fs->make<TH1F>(Form("h_E1oE7_layer%d",layer+1), Form("h_E1oE7_layer%d",layer+1), 5000, -5, 5);
	  h_E1oE19_layer[layer] = fs->make<TH1F>(Form("h_E1oE19_layer%d",layer+1), Form("h_E1oE19_layer%d",layer+1), 5000, -5, 5);
	  h_E7oE19_layer[layer] = fs->make<TH1F>(Form("h_E7oE19_layer%d",layer+1), Form("h_E7oE19_layer%d",layer+1), 5000, -5, 5);
	  h_Radius[layer] = fs->make<TH1F>(Form("h_Radius_layer%d",layer+1), Form("h_Radius_layer%d",layer+1), 5000, 0., 20.);
	}

	for(int ski=0; ski<MAXSKIROCS; ++ski) {
	  h_CM_layer[ski] = fs->make<TH1F>(Form("h_CM_skiroc%d",ski+1), "", 5000, -500, 500);
	}

	h_sum_all = fs->make<TH1F>("h_sumAll_AllLayers", "", 40010, -10, 40000);
	h_seven_all = fs->make<TH1F>("h_sum7_AllLayers", "", 40010, -10, 40000);
	h_nineteen_all = fs->make<TH1F>("h_sum19_AllLayers", "", 40010, -10, 40000);
	h_sum_all->Sumw2();
	h_seven_all->Sumw2();
	h_nineteen_all->Sumw2();
	h_sum_all_AbsW = fs->make<TH1F>("h_sumAll_AllLayers_AbsW", "", 40010, -10, 40000);
	h_seven_all_AbsW = fs->make<TH1F>("h_sum7_AllLayers_AbsW", "", 40010, -10, 40000);
	h_nineteen_all_AbsW = fs->make<TH1F>("h_sum19_AllLayers_AbsW", "", 40010, -10, 40000);
	h_sum_all_AbsW->Sumw2();
	h_seven_all_AbsW->Sumw2();
	h_nineteen_all_AbsW->Sumw2();

	HighGain_LowGain_2D = fs->make<TH2F>("h2_HGvsLG", "", 4000, 0, 4000, 4000, 0, 4000);
	Energy_LowGain_2D = fs->make<TH2F>("h2_EvsLG", "", 4000, 0, 4000, 4000, 0, 4000);

	tp_E1_vs_layer = fs->make<TProfile>("tp_E1_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E7_vs_layer = fs->make<TProfile>("tp_E7_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E19_vs_layer = fs->make<TProfile>("tp_E19_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E1oSumL_vs_layer = fs->make<TProfile>("tp_E1oSumL_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E7oSumL_vs_layer = fs->make<TProfile>("tp_E7oSumL_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E19oSumL_vs_layer = fs->make<TProfile>("tp_E19oSumL_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E1oE7_vs_layer = fs->make<TProfile>("tp_E1oE7_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E1oE19_vs_layer = fs->make<TProfile>("tp_E1oE19_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	tp_E7oE19_vs_layer = fs->make<TProfile>("tp_E7oE19_vs_layer", "", MAXLAYERS+1, 0, MAXLAYERS+1);
	

        //loading the proper weights
        if(MAXLAYERS != 8) {
	  std::cout << " update weights " << std::endl;
	  return;
	}

        for(int iL=0; iL<MAXLAYERS; ++iL){
	  if(CERN_8layers_config_ == 1){
	    Weights_8L[iL] = LayerWeight_8L_conf1[iL];
	    X0_8L[iL] = X0depth_8L_conf1[iL];
	  }
	  if(CERN_8layers_config_ == 2){
	    Weights_8L[iL] = LayerWeight_8L_conf2[iL];
	    X0_8L[iL] = X0depth_8L_conf2[iL]; 
	  }
	}


}//constructor ends here


Layer_Sum_Analyzer::~Layer_Sum_Analyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


// ------------ method called for each event  ------------
void
Layer_Sum_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

  //  std::cout << " >>> Layer_Sum_Analyzer::analyze " << std::endl;
  EVENT = (event.id()).event();

  if((event.id()).event() == 1) {
    Time_Temp = event.time().value();
    Time_Stamp[EVENT] = Time_Temp;
    Delta_Time_Stamp[EVENT] = 0.;
  } else {
    Time_Stamp[EVENT] = event.time().value();
    Delta_Time_Stamp[EVENT] = event.time().value() - Time_Temp;
    Time_Temp = event.time().value();
  }
  
  //opening Rechits
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);
  
  // looping over each rechit to fill histogram
  double commonmode[MAXSKIROCS];
  int cm_num[MAXSKIROCS];
  for(int iS=0; iS<MAXSKIROCS; ++iS){
    commonmode[iS] = 0.;
    cm_num[iS] = 0;
  }
  double max[MAXLAYERS], max_x[MAXLAYERS], max_y[MAXLAYERS];
  for(int iL=0; iL<MAXLAYERS; ++iL){
    max[iL] = max_x[iL] = max_y[iL] = 0.;
  }

  
  //CM subtraction => REPLACE WITH EXTERNAL CLASS => FIX
  for(auto Rechit : *Rechits){
    
    if(!IsCellValid.iu_iv_valid((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorsize))  continue;
   
    //getting electronics ID
    uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
    HGCalTBElectronicsId eid(EID);
    
    int n_layer = (Rechit.id()).layer() - 1;
    int n_cell_type = (Rechit.id()).cellType();
    int n_skiroc = eid.iskiroc() - 1;


    //getting X and Y coordinates
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorsize);

    HighGain_LowGain_2D->Fill(Rechit.energyLow(), Rechit.energyHigh());
    Energy_LowGain_2D->Fill(Rechit.energy(), Rechit.energyHigh());
    
    // needed? FIXME
    if(n_cell_type != 0) continue;
    
    if(Rechit.energy() > max[n_layer]) {      
      max[n_layer] = Rechit.energy();
      max_x[n_layer] = CellCentreXY.first;
      max_y[n_layer] = CellCentreXY.second;
    }

    if((Rechit.energy()) / ADCtoMIP[n_skiroc] <= CMTHRESHOLD) {
      commonmode[n_skiroc] += Rechit.energy();
      cm_num[n_skiroc]++;
    }

  }//Rechit loop ends here
  //	std::cout << " >>> found commonmode = " << commonmode << std::endl;


  for(int iS=0; iS<MAXSKIROCS; ++iS){
    commonmode[iS] = commonmode[iS]/cm_num[iS];
    //    std::cout << " value = " << commonmode[iS] << std::endl;
  }


  edm::Handle<HGCalTBRecHitCollection> Rechits1;
  event.getByToken(HGCalTBRecHitCollection_, Rechits1);

  // looping over each rechit to fill histogram
  double allcells_sum[MAXLAYERS], sevencells_sum[MAXLAYERS], nineteencells_sum[MAXLAYERS], radius[MAXLAYERS];
  double seedEnergy[MAXLAYERS];
  double x_tmp[MAXLAYERS], y_tmp[MAXLAYERS];
  int num[MAXLAYERS], sevennum[MAXLAYERS], nineteennum[MAXLAYERS];

  for(int iL=0; iL<MAXLAYERS; ++iL){
    allcells_sum[iL] = sevencells_sum[iL] = nineteencells_sum[iL] = radius[iL] = 0.;
    x_tmp[iL] = y_tmp[iL] = seedEnergy[iL] = 0.;
    num[iL] = sevennum[iL] = nineteennum[iL] = 0;
  }

  for(auto Rechit1 : *Rechits1) {

    if(!IsCellValid.iu_iv_valid((Rechit1.id()).layer(), (Rechit1.id()).sensorIU(), (Rechit1.id()).sensorIV(), (Rechit1.id()).iu(), (Rechit1.id()).iv(), sensorsize))  continue;	  
    
    //getting electronics ID
    uint32_t EID = essource_.emap_.detId2eid(Rechit1.id());
    HGCalTBElectronicsId eid(EID);
	  
    int eLayer = (Rechit1.id()).layer()-1;
    int eCellType = (Rechit1.id()).cellType();
    int eSkiroc = eid.iskiroc() - 1;


    //getting X and Y coordinates
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit1.id()).layer(), (Rechit1.id()).sensorIU(), (Rechit1.id()).sensorIV(), (Rechit1.id()).iu(), (Rechit1.id()).iv(), sensorsize);

    //FIXME >> needed?
    if(eCellType != 0) continue;

    //    std::cout << Rechit1.energy() << std::endl;

    radius[eLayer] = sqrt( pow(CellCentreXY.first - max_x[eLayer], 2) + pow(CellCentreXY.second - max_y[eLayer], 2) );

    float energyCMsub = (Rechit1.energy() - commonmode[eSkiroc]) / ADCtoMIP[eSkiroc];
    //    std::cout << "val = " << commonmode[eSkiroc] << std::endl;
    if(energyCMsub > CMTHRESHOLD){
      allcells_sum[eLayer] += energyCMsub;
      if(energyCMsub > seedEnergy[eLayer]) seedEnergy[eLayer] = energyCMsub;
    }
    
    // countin all recHits even below CM threshold
    num[eLayer]++;

    h_Radius[eLayer]->Fill(radius[eLayer]);
		
    //FIXME navigation to 7cells
    if((radius[eLayer] < maxdist && sevennum[eLayer] < 7) && (energyCMsub > CMTHRESHOLD)){
      sevencells_sum[eLayer] += energyCMsub;
      sevennum[eLayer]++;
    }

    //FIXME navigation to 19cells
    if((radius[eLayer] < 1.95 * maxdist && nineteennum[eLayer] < 19) && (energyCMsub > CMTHRESHOLD)){
//			nineteencells_sum += (LayerWeight[LAYER]*(Rechit1.energyHigh() - commonmode))/ ADCtoMIP[LAYER];
      nineteencells_sum[eLayer] += energyCMsub;
      x_tmp[eLayer] += CellCentreXY.first * energyCMsub;
      y_tmp[eLayer] += CellCentreXY.second * energyCMsub;
      nineteennum[eLayer]++;
    }
  }

  for(int iS=0; iS<MAXSKIROCS; ++iS){
    //    std::cout << "iS = " << iS << " CM = " << commonmode[iS] << std::endl; 
    h_CM_layer[iS]->Fill(commonmode[iS]);
  }
 

  float E1SumL_R = 0.;
  float E7SumL_R = 0.;
  float E19SumL_R = 0.;
  float EAllSumL_R = 0.;

  float E1SumL_AbsW = 0.;
  float E7SumL_AbsW = 0.;
  float E19SumL_AbsW = 0.;
  float EAllSumL_AbsW = 0.;

  for(int iL=0; iL<MAXLAYERS; ++iL){
    h_sum_layer[iL]->Fill(allcells_sum[iL]);
    h_Seed_layer[iL]->Fill(seedEnergy[iL]);
    h_layer_seven[iL]->Fill(sevencells_sum[iL]);
    h_layer_nineteen[iL]->Fill(nineteencells_sum[iL]);

    float layerE1 = seedEnergy[iL] * ( ADCtoMIP[iL] * weights2GeV * Weights_8L[iL] + 1.);
    float layerE7 = sevencells_sum[iL] * ( ADCtoMIP[iL] * weights2GeV * Weights_8L[iL] + 1.);
    float layerE19 = nineteencells_sum[iL] * ( ADCtoMIP[iL] * weights2GeV * Weights_8L[iL] + 1.);
    float layerEAll = allcells_sum[iL] * ( ADCtoMIP[iL] * weights2GeV * Weights_8L[iL] + 1.);

    h_sum_layer_AbsW[iL]->Fill(layerEAll);
    h_Seed_layer_AbsW[iL]->Fill(layerE1);
    h_layer_seven_AbsW[iL]->Fill(layerE7);
    h_layer_nineteen_AbsW[iL]->Fill(layerE19);
    
    tp_E1_vs_layer->Fill(iL+1, seedEnergy[iL]);
    tp_E7_vs_layer->Fill(iL+1, sevencells_sum[iL]);
    tp_E19_vs_layer->Fill(iL+1, nineteencells_sum[iL]);

    if(sevencells_sum[iL] != 0.){
      h_E1oE7_layer[iL]->Fill(seedEnergy[iL]/sevencells_sum[iL]);
      tp_E1oE7_vs_layer->Fill(iL+1, seedEnergy[iL]/sevencells_sum[iL]);
    }
    if(nineteencells_sum[iL] != 0){
      h_E1oE19_layer[iL]->Fill(seedEnergy[iL]/nineteencells_sum[iL]);
      h_E7oE19_layer[iL]->Fill(sevencells_sum[iL]/nineteencells_sum[iL]);
      tp_E1oE19_vs_layer->Fill(iL, seedEnergy[iL]/nineteencells_sum[iL]);
      tp_E7oE19_vs_layer->Fill(iL, sevencells_sum[iL]/nineteencells_sum[iL]);
    }

    h_x_layer[iL]->Fill(x_tmp[iL] / nineteencells_sum[iL]);
    h_y_layer[iL]->Fill(y_tmp[iL] / nineteencells_sum[iL]);
    h_x_y_layer[iL]->Fill(x_tmp[iL] / nineteencells_sum[iL], y_tmp[iL] / nineteencells_sum[iL]);

    E1SumL_R += seedEnergy[iL];
    E7SumL_R += sevencells_sum[iL];
    E19SumL_R += nineteencells_sum[iL];
    EAllSumL_R += allcells_sum[iL];

    E1SumL_AbsW += layerE1;
    E7SumL_AbsW += layerE7;
    E19SumL_AbsW += layerE19;
    EAllSumL_AbsW += layerEAll;
  }


  for(int iL=0; iL<MAXLAYERS; ++iL){
    if(E1SumL_R != 0.){
      tp_E1oSumL_vs_layer->Fill(iL, seedEnergy[iL]/E1SumL_R);
      tp_E7oSumL_vs_layer->Fill(iL, sevencells_sum[iL]/E7SumL_R);
      tp_E19oSumL_vs_layer->Fill(iL, nineteencells_sum[iL]/E19SumL_R);
    }
    h_sum_all->Fill(EAllSumL_R);
    h_seven_all->Fill(E7SumL_R);
    h_nineteen_all->Fill(E19SumL_R);

    h_sum_all_AbsW->Fill(EAllSumL_AbsW);
    h_seven_all_AbsW->Fill(E7SumL_AbsW);
    h_nineteen_all_AbsW->Fill(E19SumL_AbsW);
  }
  
}// analyze ends here


// ------------ method called once each job just before starting event loop  ------------
void
Layer_Sum_Analyzer::beginJob()
{
  std::cout << " Layer_Sum_Analyzer::beginJob " << std::endl;
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(mapfile_);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	}

	for(int iii = 0; iii < 16; iii++)
		ADCtoMIP[iii] = ADCtoMIP[iii] / MIP2ParticleCalib; // Converting response to 120 GeV protons to MIPs
	/*
	        for(int iii= 0; iii<16;iii++){
	            LayerWeight[iii] += 0.8;
	            LayerSumWeight += LayerWeight[iii];
	           }
	*/


}

// ------------ method called once each job just after ending the event loop  ------------
void
Layer_Sum_Analyzer::endJob()
{
  std::cout << " Layer_Sum_Analyzer::endJob " << std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Layer_Sum_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Layer_Sum_Analyzer);
