/**
	@Author: Thorben Quast
		02 June 2017
		thorben.quast@cern.ch
*/



// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
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
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
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

#include "HGCal/Reco/src/ShowerShape_2017.cc"

#define DEBUG 0


class Layer_Sum_Analyzer_2017 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit Layer_Sum_Analyzer_2017(const edm::ParameterSet&);
    ~Layer_Sum_Analyzer_2017();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    void analyze(const edm::Event& , const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    int run, event_nr;
    double beamMomentum;

    edm::EDGetToken HGCalTBRecHitCollection_;
    edm::EDGetTokenT<RunData> RunDataToken; 

    int layers_config_;

    struct {
      HGCalElectronicsMap emap_;
    } essource_;

    string mapfile_;
    int sensorsize = 128;
    
    HGCalTBTopology IsCellValid;
    HGCalTBCellVertices TheCell;



    //output quantities 
    TTree* outTree;
    double E1SumL_R;
    double E7SumL_R;
    double E19SumL_R;
    double EAllSumL_R;
    double I_Highest_R;
    double I_SecondHighest_R;
    double I_ThirdHighest_R;

    double E1SumL_R_high;
    double E7SumL_R_high;
    double E19SumL_R_high;
    double EAllSumL_R_high;
    double I_Highest_R_high;
    double I_SecondHighest_R_high;
    double I_ThirdHighest_R_high;

    double E1SumL_R_low;
    double E7SumL_R_low;
    double E19SumL_R_low;
    double EAllSumL_R_low;
    double I_Highest_R_low;
    double I_SecondHighest_R_low;
    double I_ThirdHighest_R_low;

    std::vector<double> E1SumL_R_perLayer;
    std::vector<double> E7SumL_R_perLayer;
    std::vector<double> E19SumL_R_perLayer;
    std::vector<double> EAllSumL_R_perLayer;
    std::vector<double> I_Highest_R_perLayer;
    std::vector<double> I_SecondHighest_R_perLayer;
    std::vector<double> I_ThirdHighest_R_perLayer;

    std::vector<double> E1SumL_R_high_perLayer;
    std::vector<double> E7SumL_R_high_perLayer;
    std::vector<double> E19SumL_R_high_perLayer;
    std::vector<double> EAllSumL_R_high_perLayer;
    std::vector<double> I_Highest_R_high_perLayer;
    std::vector<double> I_SecondHighest_R_high_perLayer;
    std::vector<double> I_ThirdHighest_R_high_perLayer;

    std::vector<double> E1SumL_R_low_perLayer;
    std::vector<double> E7SumL_R_low_perLayer;
    std::vector<double> E19SumL_R_low_perLayer;
    std::vector<double> EAllSumL_R_low_perLayer;
    std::vector<double> I_Highest_R_low_perLayer;
    std::vector<double> I_SecondHighest_R_low_perLayer;
    std::vector<double> I_ThirdHighest_R_low_perLayer;


    int skirocVeto;

    //fixed parameters specific for the Layer_Sum_Analyzer
    int n_layers;
    int n_skirocs;
    double* Weights_L;
    double* X0_L;
    double* ADCtoMIP;
      

    double MIP2ParticleCalib = 1;  
    double MIP2GeV_sim = 51.91e-06; //mpv muon EMM pysics list
  
    double weights2MIP = 1.;   // rescale weights from mean to MPV
    double weights2GeV = 1.e-03;

    
    //MIP indicator flags
    int nThreshold;
    double e23_over_e1;
    int e123_in_e19;

};



Layer_Sum_Analyzer_2017::Layer_Sum_Analyzer_2017(const edm::ParameterSet& iConfig) {
  
  // initialization
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  RunDataToken = consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
  layers_config_ = iConfig.getParameter<int>("layers_config");

  //configuration here
  if(layers_config_ == 1){
    n_layers = 1;
    n_skirocs = 4;
    Weights_L = new double[n_layers]; Weights_L[0] = 1.;
    X0_L = new double[n_layers]; X0_L[0] = 1.;
    ADCtoMIP = new double[n_skirocs]; ADCtoMIP[0] = 1.; ADCtoMIP[1] = 1.; ADCtoMIP[2] = 1.; ADCtoMIP[3] = 1.; 

    mapfile_ = iConfig.getParameter<std::string>("mapFile_CERN");
  } else {
    throw cms::Exception("Invalid layer configuration");
  }


  outTree = fs->make<TTree>("energySums", "energySums");
  outTree->Branch("event_nr", &event_nr, "event_nr/I");
  outTree->Branch("run", &run, "run/I");
  outTree->Branch("beamMomentum", &beamMomentum, "beamMomentum/D");
  outTree->Branch("nLayers", &n_layers, "nLayers/I");

  outTree->Branch("E1SumL_R", &E1SumL_R, "E1SumL_R/D") ;
  outTree->Branch("E7SumL_R", &E7SumL_R, "E7SumL_R/D") ;
  outTree->Branch("E19SumL_R", &E19SumL_R, "E19SumL_R/D") ;
  outTree->Branch("EAllSumL_R", &EAllSumL_R, "EAllSumL_R/D") ;
  outTree->Branch("I_Highest_R", &I_Highest_R, "I_Highest_R/D");
  outTree->Branch("I_SecondHighest_R", &I_SecondHighest_R, "I_SecondHighest_R/D");
  outTree->Branch("I_ThirdHighest_R", &I_ThirdHighest_R, "I_ThirdHighest_R/D");

  outTree->Branch("E1SumL_R_high", &E1SumL_R_high, "E1SumL_R_high/D") ;
  outTree->Branch("E7SumL_R_high", &E7SumL_R_high, "E7SumL_R_high/D") ;
  outTree->Branch("E19SumL_R_high", &E19SumL_R_high, "E19SumL_R_high/D") ;
  outTree->Branch("EAllSumL_R_high", &EAllSumL_R_high, "EAllSumL_R_high/D") ;
  outTree->Branch("I_Highest_R_high", &I_Highest_R_high, "I_Highest_R_high/D");
  outTree->Branch("I_SecondHighest_R_high", &I_SecondHighest_R_high, "I_SecondHighest_R_high/D");
  outTree->Branch("I_ThirdHighest_R_high", &I_ThirdHighest_R_high, "I_ThirdHighest_R_high/D");

  outTree->Branch("E1SumL_R_low", &E1SumL_R_low, "E1SumL_R_low/D") ;
  outTree->Branch("E7SumL_R_low", &E7SumL_R_low, "E7SumL_R_low/D") ;
  outTree->Branch("E19SumL_R_low", &E19SumL_R_low, "E19SumL_R_low/D") ;
  outTree->Branch("EAllSumL_R_low", &EAllSumL_R_low, "EAllSumL_R_low/D") ;
  outTree->Branch("I_Highest_R_low", &I_Highest_R_low, "I_Highest_R_low/D");
  outTree->Branch("I_SecondHighest_R_low", &I_SecondHighest_R_low, "I_SecondHighest_R_low/D");
  outTree->Branch("I_ThirdHighest_R_low", &I_ThirdHighest_R_low, "I_ThirdHighest_R_low/D");


  
  outTree->Branch("E1SumL_R_perLayer", &E1SumL_R_perLayer);
  outTree->Branch("E7SumL_R_perLayer", &E7SumL_R_perLayer);
  outTree->Branch("E19SumL_R_perLayer", &E19SumL_R_perLayer);
  outTree->Branch("EAllSumL_R_perLayer", &EAllSumL_R_perLayer);
  outTree->Branch("I_Highest_R_perLayer", &I_Highest_R_perLayer);
  outTree->Branch("I_SecondHighest_R_perLayer", &I_SecondHighest_R_perLayer);
  outTree->Branch("I_ThirdHighest_R_perLayer", &I_ThirdHighest_R_perLayer);

  outTree->Branch("E1SumL_R_high_perLayer", &E1SumL_R_high_perLayer);
  outTree->Branch("E7SumL_R_high_perLayer", &E7SumL_R_high_perLayer);
  outTree->Branch("E19SumL_R_high_perLayer", &E19SumL_R_high_perLayer);
  outTree->Branch("EAllSumL_R_high_perLayer", &EAllSumL_R_high_perLayer);
  outTree->Branch("I_Highest_R_high_perLayer", &I_Highest_R_high_perLayer);
  outTree->Branch("I_SecondHighest_R_high_perLayer", &I_SecondHighest_R_high_perLayer);
  outTree->Branch("I_ThirdHighest_R_high_perLayer", &I_ThirdHighest_R_high_perLayer);

  outTree->Branch("E1SumL_R_low_perLayer", &E1SumL_R_low_perLayer);
  outTree->Branch("E7SumL_R_low_perLayer", &E7SumL_R_low_perLayer);
  outTree->Branch("E19SumL_R_low_perLayer", &E19SumL_R_low_perLayer);
  outTree->Branch("EAllSumL_R_low_perLayer", &EAllSumL_R_low_perLayer);
  outTree->Branch("I_Highest_R_low_perLayer", &I_Highest_R_low_perLayer);
  outTree->Branch("I_SecondHighest_R_low_perLayer", &I_SecondHighest_R_low_perLayer);
  outTree->Branch("I_ThirdHighest_R_low_perLayer", &I_ThirdHighest_R_low_perLayer);
  


  outTree->Branch("skirocVeto", &skirocVeto, "skirocVeto/I");
  outTree->Branch("nThreshold", &nThreshold, "nThreshold/I");
  outTree->Branch("e23_over_e1", &e23_over_e1, "e23_over_e1/D");
  outTree->Branch("e123_in_e19", &e123_in_e19, "e123_in_e19/I");


}//constructor ends here


Layer_Sum_Analyzer_2017::~Layer_Sum_Analyzer_2017() {
  delete Weights_L;
  delete X0_L;
  delete ADCtoMIP;
}


// ------------ method called for each event  ------------
void Layer_Sum_Analyzer_2017::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  //if (DEBUG)  std::cout << " >>> Layer_Sum_Analyzer_2017::analyze " << std::endl;

  //opening Rechits
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);

  edm::Handle<RunData> rd;
  //get the relevant event information
  event.getByToken(RunDataToken, rd);  
  skirocVeto = (int) rd->hasDanger;
  beamMomentum = rd->energy;
  run = rd->run;
  event_nr = rd->event;


  double* max = new double[n_layers];
  double* max_x = new double[n_layers];
  double* max_y = new double[n_layers];
  
  E1SumL_R_perLayer.clear();
  E7SumL_R_perLayer.clear();
  E19SumL_R_perLayer.clear();
  EAllSumL_R_perLayer.clear();
  I_Highest_R_perLayer.clear();
  I_SecondHighest_R_perLayer.clear();
  I_ThirdHighest_R_perLayer.clear();

  E1SumL_R_high_perLayer.clear();
  E7SumL_R_high_perLayer.clear();
  E19SumL_R_high_perLayer.clear();
  EAllSumL_R_high_perLayer.clear();
  I_Highest_R_high_perLayer.clear();
  I_SecondHighest_R_high_perLayer.clear();
  I_ThirdHighest_R_high_perLayer.clear();

  E1SumL_R_low_perLayer.clear();
  E7SumL_R_low_perLayer.clear();
  E19SumL_R_low_perLayer.clear();
  EAllSumL_R_low_perLayer.clear();
  I_Highest_R_low_perLayer.clear();
  I_SecondHighest_R_low_perLayer.clear();
  I_ThirdHighest_R_low_perLayer.clear();

  for(int iL=0; iL<n_layers; ++iL){
    max[iL] = max_x[iL] = max_y[iL] = 0.;
    
    E1SumL_R_perLayer.push_back(0.);
    E7SumL_R_perLayer.push_back(0.);
    E19SumL_R_perLayer.push_back(0.);
    EAllSumL_R_perLayer.push_back(0.);
    I_Highest_R_perLayer.push_back(0.);
    I_SecondHighest_R_perLayer.push_back(0.);
    I_ThirdHighest_R_perLayer.push_back(0.);
    
    E1SumL_R_high_perLayer.push_back(0.);
    E7SumL_R_high_perLayer.push_back(0.);
    E19SumL_R_high_perLayer.push_back(0.);
    EAllSumL_R_high_perLayer.push_back(0.);
    I_Highest_R_high_perLayer.push_back(0.);
    I_SecondHighest_R_high_perLayer.push_back(0.);
    I_ThirdHighest_R_high_perLayer.push_back(0.);

    E1SumL_R_low_perLayer.push_back(0.);
    E7SumL_R_low_perLayer.push_back(0.);
    E19SumL_R_low_perLayer.push_back(0.);
    EAllSumL_R_low_perLayer.push_back(0.);
    I_Highest_R_low_perLayer.push_back(0.);
    I_SecondHighest_R_low_perLayer.push_back(0.);
    I_ThirdHighest_R_low_perLayer.push_back(0.);
  }


  for(auto Rechit : *Rechits){
    if(!IsCellValid.iu_iv_valid((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorsize))  continue;
    //getting electronics ID
    uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
    HGCalTBElectronicsId eid(EID);

    int n_layer = (Rechit.id()).layer() - 1;
    int n_cell_type = (Rechit.id()).cellType();
    int n_skiroc = (eid.iskiroc() - 1);
    if(n_cell_type != 0 && n_cell_type != 4) continue;

    double scaleFactor = 1.0;     //scale factor to account for the partial coverage by the calibration pads
    if (n_cell_type == 4) scaleFactor = 9./8.;

    //temporarily: trust low gain
    if(Rechit.energyLow() > max[n_layer]) {      
      max[n_layer] = scaleFactor * Rechit.energyLow();
      max_x[n_layer] = Rechit.getCellCenterCartesianCoordinate(0);
      max_y[n_layer] = Rechit.getCellCenterCartesianCoordinate(1);
    }


    if (DEBUG) std::cout<<"event: "<<event_nr<<"   n_layer: "<<n_layer<<"  n_cell_type: "<<n_cell_type<<"   n_skiroc: "<<n_skiroc<<std::endl;
  }//Rechit loop ends here
  
  if (DEBUG) std::cout<<"max energy: "<<max[0]<<"   max_x: "<<max_x[0]<<"  max_y: "<<max_y[0]<<std::endl;
  
  //initialisation of the shower shapes

  ShowerShape2017 shosha_tot(mapfile_, Rechits, TOT, n_layers, n_skirocs, ADCtoMIP, max_x, max_y);  
  ShowerShape2017 shosha_high(mapfile_, Rechits, ADC_HIGH, n_layers, n_skirocs, ADCtoMIP, max_x, max_y);  
  ShowerShape2017 shosha_low(mapfile_, Rechits, ADC_LOW, n_layers, n_skirocs, ADCtoMIP, max_x, max_y);  

  E1SumL_R = 0.;
  E7SumL_R = 0.;
  E19SumL_R = 0.;
  EAllSumL_R = 0.;

  E1SumL_R_high = 0.;
  E7SumL_R_high = 0.;
  E19SumL_R_high = 0.;
  EAllSumL_R_high = 0.;

  E1SumL_R_low = 0.;
  E7SumL_R_low = 0.;
  E19SumL_R_low = 0.;
  EAllSumL_R_low = 0.;


  for(int n_layer=0; n_layer<n_layers; ++n_layer){
    //baseline subtracted ADC counts
    double e1, e7, e19, eAll;
    shosha_tot.getAllEnergy(n_layer, e1, e7, e19, eAll);
    E1SumL_R += e1;
    E7SumL_R += e7;
    E19SumL_R += e19;
    EAllSumL_R += eAll;
    E1SumL_R_perLayer[n_layer] = e1;
    E7SumL_R_perLayer[n_layer] = e7;
    E19SumL_R_perLayer[n_layer] = e19;
    EAllSumL_R_perLayer[n_layer] = eAll;
    I_Highest_R_perLayer[n_layer] = shosha_tot.getCellIntensity(n_layer, 0);
    I_SecondHighest_R_perLayer[n_layer] = shosha_tot.getCellIntensity(n_layer, 1);
    I_ThirdHighest_R_perLayer[n_layer] = shosha_tot.getCellIntensity(n_layer, 2);

    //ADC counts with high gain
    double e1_high, e7_high, e19_high, eAll_high;
    shosha_high.getAllEnergy(n_layer, e1_high, e7_high, e19_high, eAll_high);
    E1SumL_R_high += e1_high;
    E7SumL_R_high += e7_high;
    E19SumL_R_high += e19_high;
    EAllSumL_R_high += eAll_high;
    E1SumL_R_high_perLayer[n_layer] = e1_high;
    E7SumL_R_high_perLayer[n_layer] = e7_high;
    E19SumL_R_high_perLayer[n_layer] = e19_high;
    EAllSumL_R_high_perLayer[n_layer] = eAll_high;
    I_Highest_R_high_perLayer[n_layer] = shosha_high.getCellIntensity(n_layer, 0);
    I_SecondHighest_R_high_perLayer[n_layer] = shosha_high.getCellIntensity(n_layer, 1);
    I_ThirdHighest_R_high_perLayer[n_layer] = shosha_high.getCellIntensity(n_layer, 2);

    //ADC counts with low gain
    double e1_low, e7_low, e19_low, eAll_low;
    shosha_low.getAllEnergy(n_layer, e1_low, e7_low, e19_low, eAll_low);
    E1SumL_R_low += e1_low;
    E7SumL_R_low += e7_low;
    E19SumL_R_low += e19_low;
    EAllSumL_R_low += eAll_low;
    E1SumL_R_low_perLayer[n_layer] = e1_low;
    E7SumL_R_low_perLayer[n_layer] = e7_low;
    E19SumL_R_low_perLayer[n_layer] = e19_low;
    EAllSumL_R_low_perLayer[n_layer] = eAll_low;
    I_Highest_R_low_perLayer[n_layer] = shosha_low.getCellIntensity(n_layer, 0);
    I_SecondHighest_R_low_perLayer[n_layer] = shosha_low.getCellIntensity(n_layer, 1);
    I_ThirdHighest_R_low_perLayer[n_layer] = shosha_low.getCellIntensity(n_layer, 2);

  }

  //assumes: only one layer to be analysed.
  I_Highest_R = shosha_tot.getCellIntensity(0, 0);
  I_SecondHighest_R = shosha_tot.getCellIntensity(0, 1);
  I_ThirdHighest_R = shosha_tot.getCellIntensity(0, 2);

  I_Highest_R_high = shosha_high.getCellIntensity(0, 0);
  I_SecondHighest_R_high = shosha_high.getCellIntensity(0, 1);
  I_ThirdHighest_R_high = shosha_high.getCellIntensity(0, 2);
  
  I_Highest_R_low = shosha_low.getCellIntensity(0, 0);
  I_SecondHighest_R_low = shosha_low.getCellIntensity(0, 1);
  I_ThirdHighest_R_low = shosha_low.getCellIntensity(0, 2);


  //implemtation of some MIP indicator flags
  double ADC_threshold = 100.;
  nThreshold = 0;
  if (I_Highest_R_low > ADC_threshold) nThreshold += 1;
  if (I_SecondHighest_R_low > ADC_threshold) nThreshold += 1;
  if (I_ThirdHighest_R_low > ADC_threshold) nThreshold += 1;
  
  e23_over_e1 = (I_SecondHighest_R_low+I_ThirdHighest_R_low)/I_Highest_R_low;
    
  e123_in_e19 = (I_Highest_R_low + I_SecondHighest_R_low + I_ThirdHighest_R_low) < E19SumL_R_low;
  
  outTree->Fill();


}// analyze ends here


// ------------ method called once each job just before starting event loop  ------------
void Layer_Sum_Analyzer_2017::beginJob() {
  if (DEBUG) std::cout << " Layer_Sum_Analyzer_2017::beginJob " << std::endl;
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
  throw cms::Exception("Unable to load electronics map");
  }

  if(layers_config_ == 1){
    for(int iskiroc = 0; iskiroc < n_skirocs; iskiroc++){
      ADCtoMIP[iskiroc] = ADCtoMIP[iskiroc] * MIP2ParticleCalib; 
      if (DEBUG) std::cout<<"ADCtoMIP["<<iskiroc<<"] = "<<ADCtoMIP[iskiroc]<<std::endl;
    }
  } else{
    throw cms::Exception("Invalid layer configuration");
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void Layer_Sum_Analyzer_2017::endJob() {
  if (DEBUG) std::cout << " Layer_Sum_Analyzer_2017::endJob " << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Layer_Sum_Analyzer_2017::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Layer_Sum_Analyzer_2017);
