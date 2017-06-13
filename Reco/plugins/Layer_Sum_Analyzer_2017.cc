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

    double Weights_L[MAXLAYERS];
    double X0_L[MAXLAYERS];
    double ADCtoMIP[MAXSKIROCS];


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

    std::vector<double> RecHits_TOT;
    std::vector<double> RecHits_HG;
    std::vector<double> RecHits_LG;

    int skirocVeto;

    //fixed parameters specific for the Layer_Sum_Analyzer
    double ADCtoMIP_CERN[MAXSKIROCS] =  {1., 1., 1., 1.};
    double MIP2ParticleCalib = 1;  
    double MIP2GeV_sim = 51.91e-06; //mpv muon EMM pysics list
  
    double weights2MIP = 1.;   // rescale weights from mean to MPV
    double weights2GeV = 1.e-03;

    double LayerWeight_1L_conf1[1] = {1.};
    double X0depth_1L_conf1[1] = {0.};
    

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

  if(layers_config_ == 1){
    for(int iL=0; iL<MAXLAYERS; ++iL){
      Weights_L[iL] = LayerWeight_1L_conf1[iL];
      X0_L[iL] = X0depth_1L_conf1[iL];
    }
    for(int iL=0; iL<MAXSKIROCS; ++iL){
      ADCtoMIP[iL] = ADCtoMIP_CERN[iL];
    }
    mapfile_ = iConfig.getParameter<std::string>("mapFile_CERN");
  } else {
    throw cms::Exception("Invalid layer configuration");
  }

  outTree = fs->make<TTree>("energySums", "energySums");
  outTree->Branch("event_nr", &event_nr, "event_nr/I");
  outTree->Branch("run", &run, "run/I");
  outTree->Branch("beamMomentum", &beamMomentum, "beamMomentum/D");

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

  outTree->Branch("RecHits_TOT", &RecHits_TOT);
  outTree->Branch("RecHits_HG", &RecHits_HG);
  outTree->Branch("RecHits_LG", &RecHits_LG);


  outTree->Branch("skirocVeto", &skirocVeto, "skirocVeto/I");
  outTree->Branch("nThreshold", &nThreshold, "nThreshold/I");
  outTree->Branch("e23_over_e1", &e23_over_e1, "e23_over_e1/D");
  outTree->Branch("e123_in_e19", &e123_in_e19, "e123_in_e19/I");


}//constructor ends here


Layer_Sum_Analyzer_2017::~Layer_Sum_Analyzer_2017() {

}


// ------------ method called for each event  ------------
void Layer_Sum_Analyzer_2017::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  if (DEBUG)  std::cout << " >>> Layer_Sum_Analyzer_2017::analyze " << std::endl;

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

  // looping over each rechit to fill histogram
  double max[MAXLAYERS], max_x[MAXLAYERS], max_y[MAXLAYERS];
  for(int iL=0; iL<MAXLAYERS; ++iL){
    max[iL] = max_x[iL] = max_y[iL] = 0.;
  }

  RecHits_TOT.clear();
  RecHits_HG.clear();
  RecHits_LG.clear();

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

    RecHits_TOT.push_back(scaleFactor * Rechit.energy());
    RecHits_HG.push_back(scaleFactor * Rechit.energyHigh());
    RecHits_LG.push_back(scaleFactor * Rechit.energyLow());

    if (DEBUG) std::cout<<"event: "<<event_nr<<"   n_layer: "<<n_layer<<"  n_cell_type: "<<n_cell_type<<"   n_skiroc: "<<n_skiroc<<std::endl;
  }//Rechit loop ends here
  
  if (DEBUG) std::cout<<"max energy: "<<max[0]<<"   max_x: "<<max_x[0]<<"  max_y: "<<max_y[0]<<std::endl;
  
  //initialisation of the shower shapes
  ShowerShape2017 shosha_tot(mapfile_, sensorsize, Rechits, TOT, ADCtoMIP, max_x, max_y);  
  ShowerShape2017 shosha_high(mapfile_, sensorsize, Rechits, ADC_HIGH, ADCtoMIP, max_x, max_y);  
  ShowerShape2017 shosha_low(mapfile_, sensorsize, Rechits, ADC_LOW, ADCtoMIP, max_x, max_y);  

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

  for(int iL=0; iL<MAXLAYERS; ++iL){
    //baseline subtracted ADC counts
    double e1, e7, e19, eAll;
    shosha_tot.getAllEnergy(iL, e1, e7, e19, eAll);
    E1SumL_R += e1;
    E7SumL_R += e7;
    E19SumL_R += e19;
    EAllSumL_R += eAll;

    //ADC counts with high gain
    double e1_high, e7_high, e19_high, eAll_high;
    shosha_high.getAllEnergy(iL, e1_high, e7_high, e19_high, eAll_high);
    E1SumL_R_high += e1_high;
    E7SumL_R_high += e7_high;
    E19SumL_R_high += e19_high;
    EAllSumL_R_high += eAll_high;

    //ADC counts with low gain
    double e1_low, e7_low, e19_low, eAll_low;
    shosha_low.getAllEnergy(iL, e1_low, e7_low, e19_low, eAll_low);
    E1SumL_R_low += e1_low;
    E7SumL_R_low += e7_low;
    E19SumL_R_low += e19_low;
    EAllSumL_R_low += eAll_low;
  }

  //assumes: only one layer to be analysed.
  I_Highest_R = shosha_tot.getCellIntensity(0);
  I_SecondHighest_R = shosha_tot.getCellIntensity(1);
  I_ThirdHighest_R = shosha_tot.getCellIntensity(2);

  I_Highest_R_high = shosha_high.getCellIntensity(0);
  I_SecondHighest_R_high = shosha_high.getCellIntensity(1);
  I_ThirdHighest_R_high = shosha_high.getCellIntensity(2);
  
  I_Highest_R_low = shosha_low.getCellIntensity(0);
  I_SecondHighest_R_low = shosha_low.getCellIntensity(1);
  I_ThirdHighest_R_low = shosha_low.getCellIntensity(2);


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
    for(int iii = 0; iii < MAXSKIROCS; iii++){
      ADCtoMIP[iii] = ADCtoMIP[iii] * MIP2ParticleCalib; // Converting response to 120 GeV protons to MIPs
      if (DEBUG) std::cout<<"ADCtoMIP["<<iii<<"] = "<<ADCtoMIP[iii]<<std::endl;
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
