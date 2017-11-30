#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "TH2Poly.h"
#include "TH1F.h"
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

#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h" //for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

#define MAXVERTICES 6
static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.

using namespace std;

class EventDisplay : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
  explicit EventDisplay(const edm::ParameterSet&);
  ~EventDisplay();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);

  // ----------member data ---------------------------
  edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;      
  edm::EDGetTokenT<RunData> RunDataToken; 
  edm::EDGetTokenT<WireChambers> DWCToken;    
  edm::EDGetTokenT<HGCalTBDWCTrack> DWCTrackToken;  
  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  int nboards;
  int sensorsize;
  std::string emapfile_;
  std::vector<int> eventsToPlot;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  edm::Service<TFileService> fs;
};

EventDisplay::EventDisplay(const edm::ParameterSet& iConfig) :
  nboards( iConfig.getUntrackedParameter<int>("NHexaBoards",17) ),
  sensorsize( iConfig.getUntrackedParameter<int>("SensorSize",128) ),
  emapfile_ (iConfig.getUntrackedParameter<std::string>("electronicsMap",""))
{
  HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
  DWCToken= consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
  DWCTrackToken= consumes<HGCalTBDWCTrack>(iConfig.getParameter<edm::InputTag>("DWCTRACKS"));

  eventsToPlot = iConfig.getParameter<std::vector<int> >("eventsToPlot");

  usesResource("TFileService");
}


EventDisplay::~EventDisplay()
{
}

void EventDisplay::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[MAXVERTICES] = {0.};
  double HexY[MAXVERTICES] = {0.};

  for(int iv = -7; iv < 8; iv++) {
    for(int iu = -7; iu < 8; iu++) {
      if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, sensorsize)) continue;
      CellXY = TheCell.GetCellCoordinatesForPlots(layerID, sensorIU, sensorIV, iu, iv, sensorsize);
      assert(CellXY.size() == 4 || CellXY.size() == 6);
      unsigned int iVertex = 0;
      for(std::vector<std::pair<double, double>>::const_iterator it = CellXY.begin(); it != CellXY.end(); it++) {
      	HexX[iVertex] =  it->first;
      	HexY[iVertex] =  it->second;
	    ++iVertex;
      }
      //Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
      poly.AddBin(CellXY.size(), HexX, HexY);
    }//loop over iu
  }//loop over iv
}
//

// ------------ method called for each event  ------------
void EventDisplay::analyze(const edm::Event& event, const edm::EventSetup& setup){
  edm::Handle<RunData> rd;
  event.getByToken(RunDataToken, rd);
 
  if (std::find(eventsToPlot.begin(), eventsToPlot.end(), rd->event)==eventsToPlot.end()) return;

  edm::Handle<HGCalTBDWCTrack> dwctrack;
  event.getByToken(DWCTrackToken, dwctrack);

  edm::Handle<WireChambers> dwcs;
  event.getByToken(DWCToken, dwcs);

  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_Token, Rechits);


  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h_Energy_board[nboards];
  TH2Poly *h_ADCHG_board[nboards];
  TH2Poly *h_ADCLG_board[nboards];
  TH2Poly *h_ADCTOT_board[nboards];
  TH2Poly *h_ADCTOA_board[nboards];

  os.str("");
  os << "Event" << rd->event;
  std::cout << "Event: " << rd->event << std::endl;
  TFileDirectory event_dir = fs->mkdir( os.str().c_str() );
  for( int iboard=0; iboard<nboards; iboard++ ){
    
    os.str("");
    os << "Energy__Board:" << iboard;
    h_Energy_board[iboard] = event_dir.make<TH2Poly>();
    h_Energy_board[iboard]->SetName( os.str().c_str() );
    h_Energy_board[iboard]->SetTitle( os.str().c_str() );
    h_Energy_board[iboard]->GetZaxis()->SetTitle( "Energy [MIP]" );
    InitTH2Poly(*h_Energy_board[iboard],iboard,0,0);

    os.str("");
    os << "HG__Board:" << iboard;
    h_ADCHG_board[iboard] = event_dir.make<TH2Poly>();
    h_ADCHG_board[iboard]->SetName( os.str().c_str() );
    h_ADCHG_board[iboard]->SetTitle( os.str().c_str() );
    h_ADCHG_board[iboard]->GetZaxis()->SetTitle( "HG [ADC]" );
    InitTH2Poly(*h_ADCHG_board[iboard],iboard,0,0);

    os.str("");
    os << "LG__Board:" << iboard;
    h_ADCLG_board[iboard] = event_dir.make<TH2Poly>();
    h_ADCLG_board[iboard]->SetName( os.str().c_str() );
    h_ADCLG_board[iboard]->SetTitle( os.str().c_str() );
    h_ADCLG_board[iboard]->GetZaxis()->SetTitle( "LG [ADC]" );
    InitTH2Poly(*h_ADCLG_board[iboard],iboard,0,0);

    os.str("");
    os << "TOT__Board:" << iboard;
    h_ADCTOT_board[iboard] = event_dir.make<TH2Poly>();
    h_ADCTOT_board[iboard]->SetName( os.str().c_str() );
    h_ADCTOT_board[iboard]->SetTitle( os.str().c_str() );
    h_ADCTOT_board[iboard]->GetZaxis()->SetTitle( "TOT [ADC]" );
    InitTH2Poly(*h_ADCTOT_board[iboard],iboard,0,0);

    os.str("");
    os << "TOA__Board:" << iboard;
    h_ADCTOA_board[iboard] = event_dir.make<TH2Poly>();
    h_ADCTOA_board[iboard]->SetName( os.str().c_str() );
    h_ADCTOA_board[iboard]->SetTitle( os.str().c_str() );
    h_ADCTOA_board[iboard]->GetZaxis()->SetTitle( "TOA [ADC]" );
    InitTH2Poly(*h_ADCTOA_board[iboard],iboard,0,0);    

  }
  std::cout<<std::endl;

  for( auto rechit : *Rechits ){
    
    HGCalTBDetId detID=rechit.id();
    if(!IsCellValid.iu_iv_valid( detID.layer(), detID.sensorIU(), detID.sensorIV(), detID.iu(), detID.iv(), sensorsize ) ) continue;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detID.layer(), detID.sensorIU(), detID.sensorIV(), detID.iu(), detID.iv(), sensorsize );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid( rechit.id().rawId() ) );
    int board = eid.iskiroc_rawhit() / 4;
    if( detID.cellType()!=1 ) {
      h_Energy_board[ board ]->Fill(iux , iuy, rechit.energy());
      h_ADCTOA_board[ board ]->Fill(iux, iuy, rechit.time());   //this is not TOA yet but the Tmax

      if (!rechit.checkFlag(HGCalTBRecHit::kHighGainSaturated))
        h_ADCHG_board[ board ]->Fill(iux, iuy, rechit.energyHigh());
      else if (!rechit.checkFlag(HGCalTBRecHit::kLowGainSaturated))
        h_ADCLG_board[ board ]->Fill(iux, iuy, rechit.energyLow());
      else if (!rechit.checkFlag(HGCalTBRecHit::kTotGainSaturated))
        h_ADCTOT_board[ board ]->Fill(iux, iuy, rechit.energyTot());
    }
  }

}

void EventDisplay::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(emapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }
}

void EventDisplay::endJob() {

}

void EventDisplay::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventDisplay);
