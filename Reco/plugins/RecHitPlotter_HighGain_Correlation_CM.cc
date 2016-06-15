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
double return_RMS_CM(double mean_sq, double mean){
      return sqrt(mean_sq - mean*mean);
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
        std::string mapfile_ = "HGCal/CondObjects/data/map_FNAL_Layer1234.txt";
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
        AllCells_Ped = fs->make<TH1F>("AllCells_Ped","AllCells_Ped",500,-250,250);
        AllCells_CM = fs->make<TH1F>("AllCells_CM","AllCells_CM",500,-250,250);
        sprintf(name, "Noise_2D_Profile_Layer");
        sprintf(title, "Noise 2D Profile Layer");
        Noise_2D_Profile = fs->make<TH2F>(name,title,512,0,511,500,-250,250);
        for(int ILayer=0;ILayer<MAXLAYERS;ILayer++){ 
	        for(int ISkiroc = 1;ISkiroc<=MAXSKIROCS;ISkiroc++){
        		for(int Channel=0; Channel<64;Channel++){
				sprintf(name, "Ski_%i_Channel_%i_Layer_%i",ISkiroc,Channel,ILayer);
                  		sprintf(title, "Ski %i Channel %i Layer %i",ISkiroc,Channel,ILayer);
            		         h_digi_layer_channel[ISkiroc-1][Channel][ILayer] = fs->make<TH1F>(name, title, 1000,-500., 500.);
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

        double Average_Pedestal_Per_Event_Full = 0;
        int Cell_counter = 0;   
        double Average_Pedestal_Per_Event_Half = 0;
        int Cell_counter_Half = 0; 
        double Average_Pedestal_Per_Event_MB = 0;
        int Cell_counter_MB = 0;
        for(auto RecHit1 : *Rechits1) {
             CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit1.id()).layer(), (RecHit1.id()).sensorIU(), (RecHit1.id()).sensorIV(), (RecHit1.id()).iu(), (RecHit1.id()).iv(), sensorsize);
//             double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
//             double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
//               Cell_counter++;
//               Average_Pedestal_Per_Event_Full += RecHit1.energyHigh();
                
               if((RecHit1.id()).cellType() == 0){              
                       Cell_counter++;
                       Average_Pedestal_Per_Event_Full += RecHit1.energyHigh();
                  }
               else if((RecHit1.id()).cellType() == 1 || (RecHit1.id()).cellType() == 2 || (RecHit1.id()).cellType() == 4){
                      Cell_counter_Half++;
                       Average_Pedestal_Per_Event_Half += RecHit1.energyHigh();
                     }
               else{
                     Cell_counter_MB++;
                     Average_Pedestal_Per_Event_MB += RecHit1.energyHigh();
                     } 

           }
	for(auto RecHit : *Rechits) {
		if(!IsCellValid.iu_iv_valid((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize))  continue;
                CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((RecHit.id()).layer(), (RecHit.id()).sensorIU(), (RecHit.id()).sensorIV(), (RecHit.id()).iu(), (RecHit.id()).iv(), sensorsize);
//                double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
//                double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
                uint32_t EID = essource_.emap_.detId2eid(RecHit.id());
                HGCalTBElectronicsId eid(EID);
                if(eid.ichan() == 23) cout<<endl<<(RecHit.id()).cellType()<<endl;
//                          TF1* Fit2= (TF1*) h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()]->GetFunction("gaus");      
                          if(!doCommonMode_CM){
                             h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()][(RecHit.id()).layer() -1]->Fill(RecHit.energyHigh());
                             Noise_2D_Profile->Fill((64*(eid.iskiroc()-1) + eid.ichan()),RecHit.energyHigh());  
                            }
                          if(doCommonMode_CM){
                             AllCells_CM->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_Full/(Cell_counter)));
                             AllCells_Ped->Fill(RecHit.energyHigh()); 
                             if((RecHit.id()).cellType() == 0) h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()][(RecHit.id()).layer() -1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_Full/(Cell_counter)));
                             else if((RecHit.id()).cellType() == 1 || (RecHit.id()).cellType() == 2 || (RecHit.id()).cellType() == 4) h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()][(RecHit.id()).layer() -1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_Half/(Cell_counter_Half)));
                             else h_digi_layer_channel[eid.iskiroc()-1][eid.ichan()][(RecHit.id()).layer() -1]->Fill(RecHit.energyHigh() - (Average_Pedestal_Per_Event_MB/(Cell_counter_MB)));
                             Noise_2D_Profile->Fill((64*(eid.iskiroc()-1) + eid.ichan()),RecHit.energyHigh() - (Average_Pedestal_Per_Event_Full/(Cell_counter)));
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
