// -*- C++ -*-
//
// Package:    HGCal/DigiPlotter_Ped_Profile
// Class:      DigiPlotter_Ped_Profile
//
/**\class DigiPlotter_Ped_Profile DigiPlotter_Ped_Profile.cc HGCal/DigiPlotter_Ped_Profile/plugins/DigiPlotter_Ped_Profile.cc

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
#include "TProfile.h"
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
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DigiPlotter_Ped_Profile : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
	explicit DigiPlotter_Ped_Profile(const edm::ParameterSet&);
	~DigiPlotter_Ped_Profile();
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
	virtual void beginJob() override;
	void analyze(const edm::Event& , const edm::EventSetup&) override;
	virtual void endJob() override;
	// ----------member data ---------------------------
	bool DEBUG = 0;
	HGCalTBTopology IsCellValid;
	HGCalTBCellVertices TheCell;
	int sensorsize = 128;// The geometry for a 256 cell sensor hasnt been implemted yet. Need a picture to do this.
	std::vector<std::pair<double, double>> CellXY;
	std::pair<double, double> CellCentreXY;
	std::vector<std::pair<double, double>>::const_iterator it;
	const static int NSAMPLES = 2;
	const static int NLAYERS  = 1;
        TH1F* Sum_Cluster_ADC;
        TH1F* CG_X;
        TH1F* CG_Y;
        TProfile* SKI1_Ped_Event;
        TProfile* SKI2_Ped_Event;       
	TH2Poly *h_digi_layer[NSAMPLES][NLAYERS][1024];
 //       TH2Poly *h_digi_layer_Raw[NSAMPLES][NLAYERS][512]; 
        TH1F* Sum_Hist_cells_SKI1[1024]; 
	const static int cellx = 15;
	const static int celly = 15;
	int Sensor_Iu = 0;
	int Sensor_Iv = 0;
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
DigiPlotter_Ped_Profile::DigiPlotter_Ped_Profile(const edm::ParameterSet& iConfig)
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
        CG_X = fs->make<TH1F>("CG_X","CG X[cm]", 100,-5,5);
        CG_Y = fs->make<TH1F>("CG_Y","CG Y[cm]", 100,-5,5);
        SKI1_Ped_Event = fs->make<TProfile>("SKI1_Ped_Event","Profile per event of pedestal for SKI1",512,1,512,0,1000);
        SKI2_Ped_Event = fs->make<TProfile>("SKI2_Ped_Event","Profile per event of pedestal for SKI2",512,1,512,0,1000);
        Sum_Cluster_ADC = fs->make<TH1F>("Sum_Cluster_ADC","Sum_Cluster_ADC",4096,0,4095); 

	for(int nsample = 0; nsample < NSAMPLES; nsample++) {
		for(int nlayers = 0; nlayers < NLAYERS; nlayers++) {
                    for(int eee = 0; eee< 1024; eee++){
//Booking a "hexagonal" histograms to display the sum of Digis for NSAMPLES, in 1 SKIROC in 1 layer. To include all layers soon. Also the 1D Digis per cell in a sensor is booked here for NSAMPLES.
			sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i", nsample, nlayers + 1,eee);
			sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i", nsample, nlayers + 1,eee);
			h_digi_layer[nsample][nlayers][eee] = fs->make<TH2Poly>();
			h_digi_layer[nsample][nlayers][eee]->SetName(name);
			h_digi_layer[nsample][nlayers][eee]->SetTitle(title);
/*
                        sprintf(name, "FullLayer_ADC%i_Layer%i_Event%i_Raw", nsample, nlayers + 1,eee);
                        sprintf(title, "Sum of adc counts per cell for ADC%i Layer%i Event%i Raw", nsample, nlayers + 1,eee);
                        h_digi_layer_Raw[nsample][nlayers][eee] = fs->make<TH2Poly>();
                        h_digi_layer_Raw[nsample][nlayers][eee]->SetName(name);
                        h_digi_layer_Raw[nsample][nlayers][eee]->SetTitle(title);
*/
                        sprintf(name, "SumCellHist_SKI1_Event%i",eee);
                        sprintf(title, "Sum of adc counts of cells  Event%i",eee);
                        Sum_Hist_cells_SKI1[eee] = fs->make<TH1F>(name,title,4096, 0 , 4095);
			for(int iv = -7; iv < 8; iv++) {
				for(int iu = -7; iu < 8; iu++) {
/*
                                        for(int eee = 0; eee< 512; eee++){
                                           }
*/
					CellXY = TheCell.GetCellCoordinatesForPlots(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
					int NumberOfCellVertices = CellXY.size();
					iii = 0;
					if(NumberOfCellVertices == 4) {
						for(it = CellXY.begin(); it != CellXY.end(); it++) {
							HalfHexX[iii] =  it->first;
							HalfHexY[iii++] =  it->second;
						}
//Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
						h_digi_layer[nsample][nlayers][eee]->AddBin(NumberOfCellVertices, HalfHexX, HalfHexY);
//                                                h_digi_layer_Raw[nsample][nlayers][eee]->AddBin(NumberOfCellVertices, HalfHexX, HalfHexY);

					} else if(NumberOfCellVertices == 6) {
						iii = 0;
						for(it = CellXY.begin(); it != CellXY.end(); it++) {
							FullHexX[iii] =  it->first;
							FullHexY[iii++] =  it->second;
						}
						h_digi_layer[nsample][nlayers][eee]->AddBin(NumberOfCellVertices, FullHexX, FullHexY);
//                                                h_digi_layer_Raw[nsample][nlayers][eee]->AddBin(NumberOfCellVertices, FullHexX, FullHexY);

					}

				}//loop over iu
			}//loop over iv
                    }//loop over number of events
		}//loop over nlayers
	}//loop over nsamples
}//contructor ends here


DigiPlotter_Ped_Profile::~DigiPlotter_Ped_Profile()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DigiPlotter_Ped_Profile::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	using namespace edm;
	using namespace std;
	std::vector<edm::Handle<SKIROC2DigiCollection> > ski;
	event.getManyByType(ski);
//        int Event = event.id().event();
	if(!ski.empty()) {

		std::vector<edm::Handle<SKIROC2DigiCollection> >::iterator i;
                double Average_Pedestal_Per_Event1_Full = 0, Average_Pedestal_Per_Event1_Half = 0, Average_Pedestal_Per_Event2_Full = 0, Average_Pedestal_Per_Event2_Half = 0;
                int Cell_counter1_Full = 0, Cell_counter2_Full = 0, Cell_counter1_Half = 0, Cell_counter2_Half = 0;
                double iux_Max=0., iyy_Max=0.;
                int ADC_TMP=0;
                int Sum_Cluster_Tmp;
                double iux_CG=0., iyy_CG=0.;
                int CG_counter=0;  
//                int counter1=0, counter2=0;
		for(i = ski.begin(); i != ski.end(); i++) {
			const SKIROC2DigiCollection& Coll = *(*i);
			if(DEBUG) cout << "SKIROC2 Digis: " << i->provenance()->branchName() << endl;
//////////////////////////////////Evaluate average pedestal per event to subtract out//////////////////////////////////
                        for(SKIROC2DigiCollection::const_iterator k = Coll.begin(); k != Coll.end(); k++) {
                                const SKIROC2DataFrame& SKI_1 = *k ;
                                int n_layer = (SKI_1.detid()).layer();
                                int n_sensor_IU = (SKI_1.detid()).sensorIU();
                                int n_sensor_IV = (SKI_1.detid()).sensorIV();
                                int n_cell_iu = (SKI_1.detid()).iu();
                                int n_cell_iv = (SKI_1.detid()).iv();
                                if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
                                if(!IsCellValid.iu_iv_valid(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize))  continue;
                                CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize);
                                double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
                                double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);
                                int nsample = 0;
                                if((SKI_1[nsample].adc() > ADC_TMP) && (SKI_1.detid()).cellType() == 0){
                                   ADC_TMP = SKI_1[nsample].adc();
                                   iux_Max = iux;
                                   iyy_Max = iyy;
                                  }
                                if((iux <= -0.25 && iux >= -3.25) && (iyy <= -0.25 && iyy>= -5.25 ) && SKI_1.detid().cellType() == 0){
                                   Cell_counter1_Full++;
                                   Average_Pedestal_Per_Event1_Full += SKI_1[nsample].adc();
                                   Sum_Hist_cells_SKI1[event.id().event() - 1]->Fill(SKI_1[nsample].adc()); 
                                  }
     
                                if((iux >= 0.25 && iux <= 3.25) && (iyy <= -0.25 && iyy>= -5.25 ) && SKI_1.detid().cellType() == 0){
                                   Cell_counter2_Full++;
                                   Average_Pedestal_Per_Event2_Full += SKI_1[nsample].adc();
                                  }

                                if((SKI_1.detid()).cellType() == 2 || (SKI_1.detid()).cellType() == 3 || (SKI_1.detid()).cellType() == 5 ){
                                    if((iux <= 0.25 && iyy>= -0.25 ) || (iux < -0.5) ){
                                          Cell_counter1_Half++;
                                          Average_Pedestal_Per_Event1_Half += SKI_1[nsample].adc();
                                        }
                                    if((iux >= 0.25 && iyy <= -0.5) || (iux>= 0.5) ){
                                         Cell_counter2_Half++;
                                         Average_Pedestal_Per_Event2_Half += SKI_1[nsample].adc();     
                                      }      
                                   }   


                        }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        SKI1_Ped_Event->Fill(event.id().event(),Average_Pedestal_Per_Event1_Full/Cell_counter1_Full);
                        SKI2_Ped_Event->Fill(event.id().event(),Average_Pedestal_Per_Event2_Full/Cell_counter2_Full);
                        iux_Max++; iyy_Max++;
                        Sum_Cluster_Tmp=0;
                        iux_CG=0;iyy_CG=0;
			for(SKIROC2DigiCollection::const_iterator j = Coll.begin(); j != Coll.end(); j++) {
                                int flag_calib = 0;
				const SKIROC2DataFrame& SKI = *j ;
				int n_layer = (SKI.detid()).layer();
				int n_sensor_IU = (SKI.detid()).sensorIU();
				int n_sensor_IV = (SKI.detid()).sensorIV();
				int n_cell_iu = (SKI.detid()).iu();
				int n_cell_iv = (SKI.detid()).iv(); 
                                if((n_cell_iu == -1 && n_cell_iv == 2) || (n_cell_iu == 2 && n_cell_iv == -4)) flag_calib =1; 
				if(DEBUG) cout << endl << " Layer = " << n_layer << " Sensor IU = " << n_sensor_IU << " Sensor IV = " << n_sensor_IV << " Cell iu = " << n_cell_iu << " Cell iu = " << n_cell_iv << endl;
				if(!IsCellValid.iu_iv_valid(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize))  continue;
				CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(n_layer, n_sensor_IU, n_sensor_IV, n_cell_iu, n_cell_iv, sensorsize);
				double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
				double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);
                                        int nsample = 0;
                                        if(flag_calib == 1){

                                            if((iux <= 0.25 && iyy>= -0.25 ) || (iux < -0.5) )
                                                h_digi_layer[nsample][n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, 0.5*(SKI[nsample].adc() - (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full))));
                                            else if((iux > -0.25 && iyy < -0.50 ) || (iux > 0.50))
                                                h_digi_layer[nsample][n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, 0.5*(SKI[nsample].adc() - (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full))));
                                          }
                                        else{
                                             if((SKI.detid()).cellType() == 2 || (SKI.detid()).cellType() == 3 || (SKI.detid()).cellType() == 5 ){
                                                   if((iux <= 0.25 && iyy>= -0.25 ) || (iux < -0.5) ){
                                                       h_digi_layer[nsample][n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (SKI[nsample].adc() - (Average_Pedestal_Per_Event1_Half/(Cell_counter1_Half))) );
                                                     }
                                                  else if((iux > -0.25 && iyy < -0.50 ) || (iux > 0.50))
                                                       h_digi_layer[nsample][n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (SKI[nsample].adc() - (Average_Pedestal_Per_Event2_Half/(Cell_counter2_Half))) );
           }
                                             else if((SKI.detid()).cellType() == 0 ){
                                                      if((iux <= 0.25 && iyy>= -0.25 ) || (iux < -0.5) ){ 
                                                          h_digi_layer[nsample][n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (SKI[nsample].adc() - (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full))) );
                                                         if(((iux < (iux_Max+2)  && iux > (iux_Max -2)) && (iyy < (iyy_Max +2) && iyy > (iyy_Max -2))) ){
//                                                            Sum_Cluster_ADC->Fill((SKI[nsample].adc() - (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full))));
                                                            Sum_Cluster_Tmp+= (SKI[nsample].adc() - (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full))); 
                                                            iux_CG += iux*(SKI[nsample].adc() - (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full)));
                                                            iyy_CG += iyy*(SKI[nsample].adc() - (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full))); 
                                                            CG_counter++;
                                                           }
                                                        }
                                                      else if((iux > -0.25 && iyy < -0.50 ) || (iux > 0.50)){
                                                               h_digi_layer[nsample][n_layer - 1][event.id().event() - 1]->Fill(iux , iyy, (SKI[nsample].adc() - (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full))) );
                                                               if((iux < (iux_Max+2)  && iux > (iux_Max -2)) && (iyy < (iyy_Max +2) && iyy > (iyy_Max -2)) ){
//                                                                   Sum_Cluster_ADC->Fill((SKI[nsample].adc() - (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full))));
                                                                   Sum_Cluster_Tmp+= (SKI[nsample].adc() - (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full)));
                                                                   iux_CG += iux*(SKI[nsample].adc() - (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full)));
                                                                   iyy_CG += iyy*(SKI[nsample].adc() - (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full)));
                                                                   CG_counter++;
                                                                 }
                                                             }
                                                     }                                  

                                            }
		}
           if((iux_Max <= 0.25 && iyy_Max>= -0.25 ) || (iux_Max < -0.5) ) Sum_Cluster_ADC->Fill(ADC_TMP - (Average_Pedestal_Per_Event1_Full/(Cell_counter1_Full)));
           else if((iux_Max > -0.25 && iyy_Max < -0.50 ) || (iux_Max > 0.50)) Sum_Cluster_ADC->Fill(ADC_TMP - (Average_Pedestal_Per_Event2_Full/(Cell_counter2_Full)));

//           Sum_Cluster_ADC->Fill(Sum_Cluster_Tmp);
           CG_X->Fill(iux_CG/Sum_Cluster_Tmp);
           CG_Y->Fill(iyy_CG/Sum_Cluster_Tmp);
         }
	} else {
		edm::LogWarning("DQM") << "No SKIROC2 Digis";
	}

}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void
DigiPlotter_Ped_Profile::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
DigiPlotter_Ped_Profile::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DigiPlotter_Ped_Profile::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiPlotter_Ped_Profile);
