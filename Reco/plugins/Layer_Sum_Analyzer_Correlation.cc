/* Need full layer, 7 cell cluster, and 19 cell cluster histograms for each layer
 * also need each for all layers summed
 * use ADC to MIP conversion of 1 MIP = 10 ADC Counts */

/**
	@Author: Ryan Quinn <ryan>
		7 July 2016
		quinn@physics.umn.edu
		
		With modifications by
                Shin-Shan Eiko Yu <syu@cern.ch> and Vieri Candelise <vieri.candelise@cern.ch>
*/



// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TF1.h"
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

double Layer_Z[16]  = {1.2, 2., 3.5, 4.3, 5.8, 6.3, 8.7, 9.5, 11.4, 12.2, 13.8, 14.6, 16.6, 17.4, 20., 20.8};
double ADCtoMIP[16] = {1.};
double LayerWeight[16] = {1.};
double LayerSumWeight = 1.;
const int CMTHRESHOLD = 30;// anything less than this value is added to the commonmode sum
const int LGCMTHRESHOLD =3; // common-mode noise for low-gain (need to be double checked, Eiko)

// applied to all layers sum after commonmode subtraction and the ADC to MIP conversion
const double PION_ALLCELLS_THRESHOLD = 15.;
const double PION_7CELLS_THRESHOLD = -100.;
const double PION_19CELLS_THRESHOLD = -100.;

const int NTYPES=6;                                           
const int NCHIPS=2;                                           
const int NCHANS=64;                                          

using namespace std;

class Layer_Sum_Analyzer_Correlation : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit Layer_Sum_Analyzer_Correlation(const edm::ParameterSet&);
	~Layer_Sum_Analyzer_Correlation();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

        bool ELECTRONS;
        bool PIONS;
 	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;

	// ----------member data ---------------------------
	edm::EDGetToken HGCalTBRecHitCollection_;
	struct {
		HGCalElectronicsMap emap_;
	} essource_;
	string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
	int sensorsize = 128;
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	HGCalTBCellVertices TheCell;
	double maxdist = (1 + sqrt (3) / 2) * HGCAL_TB_CELL::FULL_CELL_SIDE;

       // added by Eiko and Vieri                                                                                                                           
        TH2F *HighGain_LowGain_2D_lct[MAXLAYERS][NCHIPS][NTYPES];
        TProfile *pf_HighGain_LowGain_2D_cmremoved_lct[MAXLAYERS][NCHIPS][NTYPES];
        TProfile *pf_HighGain_LowGain_2D_lcc[MAXLAYERS][NCHIPS][NCHANS];
        TProfile *pf_HighGain_LowGain_2D_cmremoved_lcc[MAXLAYERS][NCHIPS][NCHANS];

	int SPILL = 0, EVENT = 0, LAYER = 0;

	map<int, double> Time_Stamp;
	map<int, double> Delta_Time_Stamp;
	double Time_Temp = 0.;
};



Layer_Sum_Analyzer_Correlation::Layer_Sum_Analyzer_Correlation(const edm::ParameterSet& iConfig)
{

  // check if it's an electron or pion beam
  if(iConfig.getParameter<std::string>("particleType")=="electron")
    {
      ELECTRONS=true;
      PIONS=false;
    }
  else if(iConfig.getParameter<std::string>("particleType")=="pion")
    {
      ELECTRONS=false;
      PIONS=true;
    }
 
  if(ELECTRONS)std::cout << "running on electron sample" << std::endl;
  if(PIONS)std::cout << "running on pion sample" << std::endl;

	// initialization
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));

	//booking the histos
	for(int layer = 0; layer < MAXLAYERS; layer++) {
		stringstream name, sevenname, nineteenname, Xname, Yname, X_Y_name;

		for(int ik = 0; ik < NCHIPS; ik++){
		  for(int ij= 0; ij < NCHANS; ij++){
		    stringstream name3;
		    name3 << "pf_HighGain_LowGain_2D_lcc" << layer+1 << Form("%02i",ik+1) << Form("%02i",ij);
		    pf_HighGain_LowGain_2D_lcc[layer][ik][ij] = fs->make<TProfile>(name3.str().c_str(), name3.str().c_str(),4000,0,4000);
		    pf_HighGain_LowGain_2D_lcc[layer][ik][ij] -> Sumw2();
  		    stringstream name4;
  		    name4 << "pf_HighGain_LowGain_2D_cmremoved_lcc" << layer+1 << Form("%02i",ik+1) << Form("%02i",ij);
  		    pf_HighGain_LowGain_2D_cmremoved_lcc[layer][ik][ij] = fs->make<TProfile>(name4.str().c_str(), name4.str().c_str(),4000,0,4000);
                    pf_HighGain_LowGain_2D_cmremoved_lcc[layer][ik][ij] -> Sumw2();
		  } // end of loop over channels, 0-63                                                                                                                                
		  for(int im= 0; im < NTYPES; im++){
		    stringstream name3;
		    name3 << "HighGain_LowGain_2D_lct" << layer+1 << Form("%02i",ik+1) << Form("%02i",im);
		    HighGain_LowGain_2D_lct[layer][ik][im] = fs->make<TH2F>(name3.str().c_str(), name3.str().c_str(),4000,0,4000,4000,0,4000);
		    HighGain_LowGain_2D_lct[layer][ik][im] -> Sumw2();

		    stringstream name4;
		    name4 << "pf_HighGain_LowGain_2D_cmremoved_lct" << layer+1 << Form("%02i",ik+1) << Form("%02i",im);
		    pf_HighGain_LowGain_2D_cmremoved_lct[layer][ik][im] = fs->make<TProfile>(name4.str().c_str(), name4.str().c_str(),4000,0,4000);
		    pf_HighGain_LowGain_2D_cmremoved_lct[layer][ik][im] -> Sumw2();

		  } // end of loop over types                                                                                                                                         
		} // end loop over skirocs   
	} // end of loop over layers
}//constructor ends here

Layer_Sum_Analyzer_Correlation::~Layer_Sum_Analyzer_Correlation()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

// ------------ method called for each event  ------------
void
Layer_Sum_Analyzer_Correlation::analyze(const edm::Event& event, const edm::EventSetup& setup)
{


	if(((event.id()).event() - 1) % (EVENTSPERSPILL * MAXLAYERS) == 0 && (event.id()).event() != 1) {
		SPILL++;
	}
	LAYER = (((event.id()).event() - 1) / EVENTSPERSPILL) % MAXLAYERS;
	EVENT = ((event.id()).event() - 1) % EVENTSPERSPILL + EVENTSPERSPILL * SPILL;
	if(EVENT == 1) {
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
	//bool FIRST(1);
	double commonmode,max, max_x, max_y;
	commonmode=max = max_x = max_y = 0.;
	int cm_num=0;

	// added by Eiko for high/low gain
	double commonmode_HG[MAXLAYERS][NCHIPS][NTYPES], commonmode_LG[MAXLAYERS][NCHIPS][NTYPES];
	int cm_num_HG[MAXLAYERS][NCHIPS][NTYPES];
        int cm_num_LG[MAXLAYERS][NCHIPS][NTYPES];
	
	for(int il=0; il<MAXLAYERS; il++){
	  for(int ic=0; ic<NCHIPS; ic++){
	    for(int it=0; it<NTYPES; it++){
	      commonmode_HG[il][ic][it]=0.;
	      commonmode_LG[il][ic][it]=0.;
	      cm_num_HG[il][ic][it]=0;
	      cm_num_LG[il][ic][it]=0;
	    }	    
	  }
	}

	for(auto Rechit : *Rechits) {

		//getting electronics ID
		uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
		HGCalTBElectronicsId eid(EID);

		//getting X and Y coordinates
		CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), sensorsize);
		int type = (Rechit.id()).cellType();
		int n_layer = (Rechit.id()).layer()-1;
		int skiroc_chip = (eid.iskiroc()-1)%2;
                int chan = eid.ichan();
		//std::cout << "n_layer = " << n_layer << "\t" << type << "\t" << skiroc_chip << "\t" << chan << std::endl;                                                    
  
                HighGain_LowGain_2D_lct[n_layer][skiroc_chip][type]->Fill(Rechit.energyLow(),Rechit.energyHigh());
                pf_HighGain_LowGain_2D_lcc[n_layer][skiroc_chip][chan]->Fill(Rechit.energyLow(),Rechit.energyHigh());
 
		if((Rechit.id()).cellType() != 0) continue;

		if((Rechit.energyHigh()) / ADCtoMIP[LAYER] <= CMTHRESHOLD) {		  
		  commonmode_HG[n_layer][skiroc_chip][type] += Rechit.energyHigh();
		  cm_num_HG[n_layer][skiroc_chip][type]++;
		  commonmode += Rechit.energyHigh();
		  cm_num++;
		}
		if((Rechit.energyLow())/ADCtoMIP[LAYER] <= LGCMTHRESHOLD){
		  commonmode_LG[n_layer][skiroc_chip][type] += Rechit.energyLow();
		  cm_num_LG[n_layer][skiroc_chip][type]++;
                }
	}//Rechit loop ends here

	commonmode = cm_num ==0? 0: commonmode/cm_num;
	for(int il=0; il<MAXLAYERS; il++){
	  for(int ic=0; ic<NCHIPS; ic++){
	    for(int it=0; it<NTYPES; it++){
	      commonmode_HG[il][ic][it] = cm_num_HG[il][ic][it]==0? 0: commonmode_HG[il][ic][it]/cm_num_HG[il][ic][it];
	      commonmode_LG[il][ic][it] = cm_num_LG[il][ic][it]==0? 0: commonmode_LG[il][ic][it]/cm_num_LG[il][ic][it]; // added by Eiko	
	      // std::cout << "common mode_HG[" << il << "]["<< ic << "][" << it << "] = " << commonmode_HG[il][ic][it] << std::endl; 
	      // std::cout << "common mode_LG[" << il << "][" << ic << "][" << it << "] = " << commonmode_LG[il][ic][it] << std::endl;                                               
	    }
	  }
	}

	//      now plot histograms after subtracting common mode                                                                                                               
	for(auto Rechit : *Rechits){
	  uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
	  HGCalTBElectronicsId eid(EID);
	  int type = (Rechit.id()).cellType();
	  int n_layer = (Rechit.id()).layer()-1;
	  int skiroc_chip = (eid.iskiroc()-1)%2;
	  int chan = eid.ichan();
	  pf_HighGain_LowGain_2D_cmremoved_lct[n_layer][skiroc_chip][type]->Fill(Rechit.energyLow()-commonmode_LG[n_layer][skiroc_chip][type],Rechit.energyHigh()-commonmode_HG[n_layer][skiroc_chip][type]);
	  pf_HighGain_LowGain_2D_cmremoved_lcc[n_layer][skiroc_chip][chan]->Fill(Rechit.energyLow()-commonmode_LG[n_layer][skiroc_chip][type], Rechit.energyHigh()-commonmode_HG[n_layer][skiroc_chip][type]);
	}
	
	edm::Handle<HGCalTBRecHitCollection> Rechits1;
	event.getByToken(HGCalTBRecHitCollection_, Rechits1);

	// looping over each rechit to fill histogram
	//double x_tmp = 0., y_tmp = 0.;
	//int num, sevennum, nineteennum;
	//num = sevennum = nineteennum = 0;

}// analyze ends here


// ------------ method called once each job just before starting event loop  ------------
void
Layer_Sum_Analyzer_Correlation::beginJob()
{
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(mapfile_);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	}

	for(int iii = 0; iii < 16; iii++)
		ADCtoMIP[iii] = ADCtoMIP[iii] / 1.3; // Converting response to 120 GeV protons to MIPs
}

// ------------ method called once each job just after ending the event loop  ------------
void
Layer_Sum_Analyzer_Correlation::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Layer_Sum_Analyzer_Correlation::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Layer_Sum_Analyzer_Correlation);
