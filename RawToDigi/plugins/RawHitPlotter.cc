#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TTree.h"
#include <fstream>
#include <sstream>
#include <algorithm>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/DataFormats/interface/HGCalTBRawHitCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h" 
#include <iomanip>
#include <set>

#include "HGCal/Reco/interface/PulseFitter.h"
double parabolicFit(std::vector<double> x, std::vector<double> y) {
      if (x.size()!=3) return -1;

      //energy reconstruction from parabolic function a*x^2 + b*x + c
      double _a = ( (y[2]-y[1])/(x[2]-x[1]) - (y[1]-y[0])/(x[1]-x[0]) )/( (pow(x[2],2)-pow(x[1],2))/(x[2]-x[1]) - (pow(x[1],2)-pow(x[0],2))/(x[1]-x[0]) );
      double _b = (y[2]-y[1]+_a*(pow(x[1],2)-pow(x[2],2)))/(x[2]-x[1]);
      //double _c = y[0]-_a*pow(x[0],2)-_b*x[0];

      double max_x = (_a < 0) ? -_b/(2*_a) : 0;   //require maximum <--> a<0, unit is ns        
      double f_x = (max_x<=150. && max_x>=25.) ? _a*pow(max_x,2)+_b*max_x : -1000.;

      return f_x;
}


#define MAXVERTICES 6
static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.

struct channelInfo{
  channelInfo(){;}
  void init(){
    for(int i=0; i<NUMBER_OF_TIME_SAMPLES; i++){
      meanHGMap[i]=0;
      meanLGMap[i]=0;
      rmsHGMap[i]=0;
      rmsLGMap[i]=0;
      counterHGMap[i]=0;
      counterLGMap[i]=0;
    }
  }
  int key;
  float meanHGMap[NUMBER_OF_TIME_SAMPLES];
  float meanLGMap[NUMBER_OF_TIME_SAMPLES];
  float rmsHGMap[NUMBER_OF_TIME_SAMPLES];
  float rmsLGMap[NUMBER_OF_TIME_SAMPLES];
  int counterHGMap[NUMBER_OF_TIME_SAMPLES];
  int counterLGMap[NUMBER_OF_TIME_SAMPLES];
  TH1F* h_adcHigh[NUMBER_OF_TIME_SAMPLES];
  TH1F* h_adcLow[NUMBER_OF_TIME_SAMPLES];
};

class RawHitPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit RawHitPlotter(const edm::ParameterSet&);
  ~RawHitPlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);

  std::string m_electronicMap;

  edm::EDGetTokenT<WireChambers> MWCToken;
  edm::EDGetTokenT<RunData> RunDataToken; 

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  int m_sensorsize;
  bool m_eventPlotter;
  bool m_subtractCommonMode;
  double m_commonModeThreshold; //currently not use (need to implement the "average" option in CommonMode.cc)

  int m_evtID;
  uint16_t m_numberOfBoards;

  std::map<int,channelInfo*> m_channelMap;

  std::map<int,TH1F*> m_h_cmHigh;
  std::map<int,TH1F*> m_h_cmLow;

  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;
  
  TTree* recHitsTree;
  int tree_run;
  int tree_board;
  int tree_skiroc;
  int tree_channel;
  double tree_lg3, tree_lg0;
  double tree_hg3, tree_hg0;
  double tree_lg3_cms, tree_lg0_cms;
  double tree_hg3_cms, tree_hg0_cms; 
  double lowGain_fit, highGain_fit, lowGain_cm_fit, highGain_cm_fit;
  double tree_hg3_cm, tree_hg0_cm, tree_lg3_cm, tree_lg0_cm;
  double tree_dwc_x_at_0, tree_dwc_y_at_0;
  double tree_dwc_mx, tree_dwc_bx, tree_dwc_my, tree_dwc_by;

};

RawHitPlotter::RawHitPlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_subtractCommonMode(iConfig.getUntrackedParameter<bool>("SubtractCommonMode",false))
  {
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));
  MWCToken = consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
  RunDataToken = consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));


  m_evtID=0;
  
  std::cout << iConfig.dump() << std::endl;
}


RawHitPlotter::~RawHitPlotter()
{

}

void RawHitPlotter::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  std::ostringstream os( std::ostringstream::ate );
  TH1F* htmp1;
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){ 
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      TFileDirectory adcdir = dir.mkdir( "ADC" );
      int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*(iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iski)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
      for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
  HGCalTBElectronicsId eid(skiId,ichan);      
  if( !essource_.emap_.existsEId(eid) ) continue;
  int key=ib*1000+iski*100+ichan;
  channelInfo *cif=new channelInfo();
  cif->key=key;
  cif->init();
  for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
    os.str("");
    os << "HighGain_Channel" << ichan << "_TS" << it ;
    htmp1=adcdir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
    cif->h_adcHigh[it]=htmp1;

    os.str("");
    os << "LowGain_Channel" << ichan << "_TS" << it ;
    htmp1=adcdir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
    cif->h_adcLow[it]=htmp1;
  }
  m_channelMap.insert( std::pair<int,channelInfo*>(key,cif) );
      }
      TFileDirectory cmdir = dir.mkdir( "CommonMode" );
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
  os.str("");
  os << "HighGain_TS" << it ;
  htmp1=cmdir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
  m_h_cmHigh.insert( std::pair<int,TH1F*>(ib*1000+iski*100+it, htmp1) );

  os.str("");
  os << "LowGain_TS" << it ;
  htmp1=cmdir.make<TH1F>(os.str().c_str(),os.str().c_str(),4000,-500,3500);
  m_h_cmLow.insert( std::pair<int,TH1F*>(ib*1000+iski*100+it, htmp1) );
      }
    }
  }

  recHitsTree = fs->make<TTree>("rawHitsTree", "rawHitsTree");
  recHitsTree->Branch("run", &tree_run);
  recHitsTree->Branch("board", &tree_board);
  recHitsTree->Branch("skiroc", &tree_skiroc);
  recHitsTree->Branch("channel", &tree_channel);
  recHitsTree->Branch("lg3", &tree_lg3);
  recHitsTree->Branch("hg3", &tree_hg3);
  recHitsTree->Branch("lg0", &tree_lg0);
  recHitsTree->Branch("hg0", &tree_hg0);
  recHitsTree->Branch("lg3_cms", &tree_lg3_cms);
  recHitsTree->Branch("hg3_cms", &tree_hg3_cms);
  recHitsTree->Branch("lg0_cms", &tree_lg0_cms);
  recHitsTree->Branch("hg0_cms", &tree_hg0_cms);
  recHitsTree->Branch("lg3_cm", &tree_lg3_cm);
  recHitsTree->Branch("hg3_cm", &tree_hg3_cm);
  recHitsTree->Branch("lg0_cm", &tree_lg0_cm);
  recHitsTree->Branch("hg0_cm", &tree_hg0_cm);
  recHitsTree->Branch("lowGain_fit", &lowGain_fit);
  recHitsTree->Branch("highGain_fit", &highGain_fit);
  recHitsTree->Branch("lowGain_cm_fit", &lowGain_cm_fit);
  recHitsTree->Branch("highGain_cm_fit", &highGain_cm_fit);
  recHitsTree->Branch("dwc_x_at_0", &tree_dwc_x_at_0);
  recHitsTree->Branch("dwc_y_at_0", &tree_dwc_y_at_0);
  recHitsTree->Branch("dwc_mx", &tree_dwc_mx);
  recHitsTree->Branch("dwc_bx", &tree_dwc_bx);
  recHitsTree->Branch("dwc_my", &tree_dwc_my);
  recHitsTree->Branch("dwc_by", &tree_dwc_by);
}

void RawHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<RunData> rd;
  event.getByToken(RunDataToken, rd);
  tree_run = rd->run;

  //get the multi wire chambers
  edm::Handle<WireChambers> dwcs;
  event.getByToken(MWCToken, dwcs);

  tree_dwc_x_at_0 = tree_dwc_y_at_0 = -999.;
  tree_dwc_mx = tree_dwc_my = -999.;
  tree_dwc_bx = tree_dwc_by = -999.;
  
  if (dwcs->at(0).goodMeasurement) {
    tree_dwc_mx = tree_dwc_my = 0.;
    tree_dwc_bx = dwcs->at(0).x;
    tree_dwc_by = dwcs->at(0).y;
  } else if (dwcs->at(1).goodMeasurement) {
    tree_dwc_mx = tree_dwc_my = 0.;
    tree_dwc_bx = dwcs->at(1).x;
    tree_dwc_by = dwcs->at(1).y;
  }

  if (dwcs->at(0).goodMeasurement && dwcs->at(1).goodMeasurement){
    double x1 = dwcs->at(0).x;
    double x2 = dwcs->at(1).x;
    double y1 = dwcs->at(0).y;
    double y2 = dwcs->at(1).y;
    double z1 = dwcs->at(0).z;
    double z2 = dwcs->at(1).z;    

    tree_dwc_mx = (x2-x1)/(z2-z1);
    tree_dwc_my = (y2-y1)/(z2-z1);
    tree_dwc_bx = x1 - tree_dwc_mx * z1 ; 
    tree_dwc_by = y1 - tree_dwc_my * z1 ; 
  } 

  tree_dwc_x_at_0 = tree_dwc_bx;
  tree_dwc_y_at_0 = tree_dwc_by;



  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  if( !hits->size() )
    return;

  std::map<int,TH2Poly*>  polyMap;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
  TH2Poly *h=dir.make<TH2Poly>();
  os.str("");
  os<<"HexaBoard"<<ib<<"_TimeSample"<<it;
  h->SetName(os.str().c_str());
  h->SetTitle(os.str().c_str());
  InitTH2Poly(*h, (int)ib, 0, 0);
  polyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      }
    }
  }
  
  CommonMode cm(essource_.emap_); //default is common mode per chip using the median
  cm.Evaluate( hits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();
  

  for( std::map<int,commonModeNoise>::iterator it=cmMap.begin(); it!=cmMap.end(); ++it ){
    int iski=(it->first-1)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ilayer=(it->first-1)/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int key=ilayer*1000+iski*100;
    for( uint16_t ts=0; ts<NUMBER_OF_TIME_SAMPLES; ts++ ){
      //m_h_cmHigh[key]->Fill( it->second.fullHG[ts] );
      //m_h_cmLow[key]->Fill( it->second.fullLG[ts] );
      key+=1;
    }
  }
  
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iboard=hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ichan=hit.channel();
    int iski=hit.skiroc();
    std::pair<int,HGCalTBDetId> p( iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan,hit.detid() );
    setOfConnectedDetId.insert(p);
    channelInfo* cif=m_channelMap[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan];
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float highGain,lowGain;
      if( m_subtractCommonMode ){
    iski = eid.iskiroc();
    float subHG(0),subLG(0);
    switch ( hit.detid().cellType() ){
    case 0 : subHG=cmMap[iski].fullHG[it]; subLG=cmMap[iski].fullLG[it]; break;
    case 2 : subHG=cmMap[iski].halfHG[it]; subLG=cmMap[iski].halfLG[it]; break;
    case 3 : subHG=cmMap[iski].mouseBiteHG[it]; subLG=cmMap[iski].mouseBiteLG[it]; break;
    case 4 : subHG=cmMap[iski].outerHG[it]; subLG=cmMap[iski].outerLG[it]; break;
    }
    highGain=hit.highGainADC(it)-subHG;
    lowGain=hit.lowGainADC(it)-subLG;
      }
      else{
    highGain=hit.highGainADC(it);
    lowGain=hit.lowGainADC(it);
      }
      iski=hit.skiroc();
      if( !hit.isUnderSaturationForHighGain() ){
  cif->h_adcHigh[it]->Fill(highGain);
  cif->meanHGMap[it]+=highGain;
  cif->rmsHGMap[it]+=highGain*highGain;
  cif->counterHGMap[it]+=1;
      }
      if( !hit.isUnderSaturationForLowGain() ){
  cif->h_adcLow[it]->Fill(lowGain);
  cif->meanLGMap[it]+=lowGain;
  cif->rmsLGMap[it]+=lowGain*lowGain;
  cif->counterLGMap[it]+=1;
      }
      if(!m_eventPlotter||!IsCellValid.iu_iv_valid(hit.detid().layer(),hit.detid().sensorIU(),hit.detid().sensorIV(),hit.detid().iu(),hit.detid().iv(),m_sensorsize))  continue;
      CellCentreXY=TheCell.GetCellCentreCoordinatesForPlots(hit.detid().layer(),hit.detid().sensorIU(),hit.detid().sensorIV(),hit.detid().iu(),hit.detid().iv(),m_sensorsize);
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
      polyMap[ 100*iboard+it ]->Fill(iux/2 , iuy, highGain);
    }

    if (hit.detid().cellType()==0)  {  //only full cells for now
      std::vector<double> sampleLG, sampleHG, sampleLGCM, sampleHGCM, sampleT;

      
      //seven time samples for fitting
      for (int i=0; i<7; i++) {
        sampleT.push_back(i*25.);
        sampleLG.push_back(hit.lowGainADC(i));
        sampleHG.push_back(hit.highGainADC(i));
        sampleLGCM.push_back(hit.lowGainADC(i)-cmMap[tree_skiroc].fullHG[i]);
        sampleHGCM.push_back(hit.highGainADC(i)-cmMap[tree_skiroc].fullHG[i]);
      }

      tree_board = hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
      tree_skiroc = hit.skiroc();
      tree_channel = hit.channel();
      
      tree_lg3 = hit.lowGainADC(3);
      tree_hg3 = hit.highGainADC(3);
      tree_lg3_cms = hit.lowGainADC(3)-cmMap[tree_skiroc].fullLG[3];
      tree_hg3_cms = hit.highGainADC(3)-cmMap[tree_skiroc].fullHG[3];

      tree_lg0 = hit.lowGainADC(0);
      tree_hg0 = hit.highGainADC(0);
      tree_lg0_cms = hit.lowGainADC(0)-cmMap[tree_skiroc].fullLG[0];
      tree_hg0_cms = hit.highGainADC(0)-cmMap[tree_skiroc].fullHG[0];
    
      //pulse fitting
      PulseFitter fitter(0,150);
      lowGain_cm_fit = highGain_cm_fit = lowGain_fit = highGain_fit = 9999.;
      
      //high gain fit
      float en0=hit.highGainADC(0)-cmMap[tree_skiroc].fullHG[0];
      float en3=hit.highGainADC(3)-cmMap[tree_skiroc].fullHG[3];
      float en6=hit.highGainADC(6)-cmMap[tree_skiroc].fullHG[6];
      if( en0<en3 && en3>en6){
        PulseFitterResult fithgcm;
        fitter.run(sampleT, sampleHGCM, fithgcm);
        if (fithgcm.status==0) {
          highGain_cm_fit = fithgcm.amplitude;
        }
      }

      //low gain fit
      en0=hit.lowGainADC(0)-cmMap[tree_skiroc].fullLG[0];
      en3=hit.lowGainADC(3)-cmMap[tree_skiroc].fullLG[3];
      en6=hit.lowGainADC(6)-cmMap[tree_skiroc].fullLG[6];
      if( en0<en3 && en3>en6){
        PulseFitterResult fitlgcm;
        fitter.run(sampleT, sampleLGCM, fitlgcm);
        if (fitlgcm.status==0) {
          lowGain_cm_fit = fitlgcm.amplitude;
        }
      }

      tree_lg3_cm = cmMap[tree_skiroc].fullLG[3];
      tree_hg3_cm = cmMap[tree_skiroc].fullHG[3];
      tree_lg0_cm = cmMap[tree_skiroc].fullLG[0];
      tree_hg0_cm = cmMap[tree_skiroc].fullHG[0];

      recHitsTree->Fill();
    }
  }
}

void RawHitPlotter::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[MAXVERTICES] = {0.};
  double HexY[MAXVERTICES] = {0.};
  for(int iv = -7; iv < 8; iv++) {
    for(int iu = -7; iu < 8; iu++) {
      if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize)) continue;
      CellXY = TheCell.GetCellCoordinatesForPlots(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize);
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


void RawHitPlotter::endJob()
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  std::map<int,TH2Poly*>  hgMeanMap;
  std::map<int,TH2Poly*>  lgMeanMap;
  std::map<int,TH2Poly*>  hgRMSMap;
  std::map<int,TH2Poly*>  lgRMSMap;
  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h;
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    os.str("");
    os << "HexaBoard" << ib ;
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    TFileDirectory hgpdir = dir.mkdir( "HighGainPedestal" );
    TFileDirectory lgpdir = dir.mkdir( "LowGainPedestal" );
    TFileDirectory hgndir = dir.mkdir( "HighGainNoise" );
    TFileDirectory lgndir = dir.mkdir( "LowGainNoise" );
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      h=hgpdir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      hgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgpdir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      lgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=hgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      hgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      lgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
    }
  }
  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA) ;
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
    channelInfo* cif=m_channelMap[iboard*1000+iski*100+ichan];
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float hgMean=cif->meanHGMap[it]/cif->counterHGMap[it];
      float lgMean=cif->meanLGMap[it]/cif->counterLGMap[it];
      float hgRMS=std::sqrt(cif->rmsHGMap[it]/cif->counterHGMap[it]-cif->meanHGMap[it]/cif->counterHGMap[it]*cif->meanHGMap[it]/cif->counterHGMap[it]);
      float lgRMS=std::sqrt(cif->rmsLGMap[it]/cif->counterLGMap[it]-cif->meanLGMap[it]/cif->counterLGMap[it]*cif->meanLGMap[it]/cif->counterLGMap[it]);
      hgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgMean );
      lgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgMean );
      hgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgRMS );
      lgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgRMS );
    }
  }
}

void RawHitPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RawHitPlotter);