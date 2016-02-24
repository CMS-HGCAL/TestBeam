// -*- C++ -*-
//
// Package:    HGCal/RecHitPlotter
// Class:      RecHitPlotter
// 
/**\class RecHitPlotter RecHitPlotter.cc HGCal/RecHitPlotter/plugins/RecHitPlotter.cc

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
using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class RecHitPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit RecHitPlotter(const edm::ParameterSet&);
      ~RecHitPlotter();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
              void analyze(const edm::Event& , const edm::EventSetup&) override; 
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetToken HGCalTBRecHitCollection_;
      HGCalTBTopology IsCellValid;
      HGCalTBCellVertices TheCell;
     int sensorsize = 128;
     std::vector<std::pair<double,double>> CellXY;
     std::pair<double,double> CellCentreXY;
     std::vector<std::pair<double,double>>::const_iterator it;

      TH2Poly *h_digi_layer,*h_digi_layer_Occupancy;
      const static int cellx = 15;
      const static int celly = 15;   
      TH1F  *h_digi_layer_cell[cellx][celly];
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
RecHitPlotter::RecHitPlotter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   sprintf(name,"FullLayer");
   sprintf(title,"Sum of RecHit Energy");
   h_digi_layer = fs->make<TH2Poly>();
   h_digi_layer->SetName(name);
   h_digi_layer->SetTitle(title);
   sprintf(name,"Occupancy");
   sprintf(title,"Occupancy");
   h_digi_layer_Occupancy = fs->make<TH2Poly>();
   h_digi_layer_Occupancy->SetName(name);
   h_digi_layer_Occupancy->SetTitle(title);
   HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));

   const int HalfHexVertices=4;
   double HalfHexX[HalfHexVertices]={0.};
   double HalfHexY[HalfHexVertices]={0.};
   const int FullHexVertices=6;
   double FullHexX[FullHexVertices]={0.};
   double FullHexY[FullHexVertices]={0.};
   int iii=0;
   for(int iv=-7;iv<8;iv++){
      for(int ix=-7;ix<8;ix++){
         if(!IsCellValid.ix_iv_valid(ix,iv,sensorsize)) continue;
         sprintf(name,"Cell_X_%i_V_%i",ix,iv);
         sprintf(title,"Rechits for Cell_X_%i_V_%i",ix,iv);
         h_digi_layer_cell[ix+7][iv+7] = fs->make<TH1F>(name,title,100,0.,40.);
         h_digi_layer_cell[ix+7][iv+7]->GetXaxis()->SetTitle("RecHits[GeV]");        
         CellXY = TheCell.GetCellCoordinates(ix,iv, sensorsize);
         int NumberOfCellVertices = CellXY.size();
         iii = 0;
         if(NumberOfCellVertices == 4){
            for(it=CellXY.begin();it != CellXY.end();it++){
                HalfHexX[iii] =  it->first;
                HalfHexY[iii++] =  it->second;
               }
            h_digi_layer->AddBin(NumberOfCellVertices,HalfHexX, HalfHexY);
            h_digi_layer_Occupancy->AddBin(NumberOfCellVertices,HalfHexX, HalfHexY);
           }
         else if(NumberOfCellVertices == 6){
            iii=0;
            for(it=CellXY.begin();it != CellXY.end();it++){
                FullHexX[iii] =  it->first;
                FullHexY[iii++] =  it->second;
               }
            h_digi_layer->AddBin(NumberOfCellVertices,FullHexX, FullHexY);
            h_digi_layer_Occupancy->AddBin(NumberOfCellVertices,FullHexX, FullHexY);
           }

        }
      }
}//contructor ends here


RecHitPlotter::~RecHitPlotter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RecHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;
   using namespace std;

   edm::Handle<HGCalTBRecHitCollection> Rechits;
   event.getByToken(HGCalTBRecHitCollection_, Rechits);
  
   for(auto RecHit: *Rechits){
//       CellXY = TheCell.GetCellCoordinates((RecHit.id()).ix(),(RecHit.id()).iv(), sensorsize);
       CellCentreXY = TheCell.GetCellCentreCoordinates((RecHit.id()).ix(),(RecHit.id()).iv(), sensorsize);
       double ixx = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
       double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);
       if(CellCentreXY.first != -1){
           h_digi_layer->Fill(ixx , iyy,RecHit.energy());
           if(RecHit.energy() > 5) h_digi_layer_Occupancy->Fill(ixx , iyy, 0);
         }
        
       if(IsCellValid.ix_iv_valid((RecHit.id()).ix(),(RecHit.id()).iv(),sensorsize)) h_digi_layer_cell[7 + (RecHit.id()).ix()][7 + (RecHit.id()).iv()]->Fill(RecHit.energy());        
//       if(CellXY.size() > 0) std::cout<<std::endl<<"Vertex 1 x = "<< (CellXY.begin())->first;
     }
  
 

}//analyze method ends here


// ------------ method called once each job just before starting event loop  ------------
void 
RecHitPlotter::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecHitPlotter::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitPlotter);
