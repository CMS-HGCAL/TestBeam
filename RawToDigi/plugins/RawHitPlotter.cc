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
#include <iomanip>
#include <set>

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

  struct {
    HGCalElectronicsMap emap_;
  } essource_;

  int m_NHexaBoards;

  int m_sensorsize;
  bool m_eventPlotter;
  bool m_subtractCommonMode;
  double m_commonModeThreshold; //currently not use (need to implement the "average" option in CommonMode.cc)

  int m_evtID;
  uint16_t m_numberOfBoards;

  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;

  std::map<int,channelInfo*> m_channelMap;

  std::map<int, TH2F*> m_h_HighVsLowGainTS3;
  std::map<int, TH2F*> m_h_LowGainVsTOTTS3;

  std::map<int,TH1F*> m_h_cmHigh;
  std::map<int,TH1F*> m_h_cmLow;
  std::map<int,TH2F*> m_h_cmHighTime;
  std::map<int,TH2F*> m_h_cmLowTime;
  std::map<int,TH1F*> m_h_nhitPerSkirocUnderSat;

  TH1F* m_h_nhitUnderSat;
  TH2F* m_h_nhitUnderSat_VS_SkirocID;
  TH1F* m_h_EEQuality;
  TH1F* m_h_FHQuality;
  TH1F* m_h_CenterQuality;
  TH1F* m_h_FullQuality;
  TTree* m_tree;
  int m_layerID;
  int m_skirocID;
  float m_cmHigh[NUMBER_OF_TIME_SAMPLES];
  float m_cmLow[NUMBER_OF_TIME_SAMPLES];
  int m_time[NUMBER_OF_TIME_SAMPLES];
};

RawHitPlotter::RawHitPlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_NHexaBoards(iConfig.getUntrackedParameter<int>("NHexaBoards", 10)), 
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_subtractCommonMode(iConfig.getUntrackedParameter<bool>("SubtractCommonMode",false))
  {
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));
  
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
  TH2F* htmp2;
  for(int ib = 0; ib<m_NHexaBoards; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){ 
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      TFileDirectory adcdir = dir.mkdir( "ADC" );
      int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*ib+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iski)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
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
      htmp2=cmdir.make<TH2F>("HighGain_RatiosVsTS","HighGain_RatiosVsTS",12,-0.5,11.5,1000,-10,10);
      m_h_cmHighTime.insert( std::pair<int,TH2F*>(ib*10+iski, htmp2) );
      htmp2=cmdir.make<TH2F>("LowGain_RatiosVsTS","LowGain_RatiosVsTS",12,-0.5,11.5,1000,-10,10);
      m_h_cmLowTime.insert( std::pair<int,TH2F*>(ib*10+iski, htmp2) );
      htmp1=cmdir.make<TH1F>("NhitUnderZero","NhitUnderZero",33,0,33);
      m_h_nhitPerSkirocUnderSat.insert( std::pair<int,TH1F*>(ib*10+iski, htmp1) );
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
    
      TFileDirectory gaindir = dir.mkdir( "Gains" );
      htmp2=gaindir.make<TH2F>("HighGainVsLowGainTS3","HighGainVsLowGainTS3",200,-500.,1500,200,-500.,3500);
      m_h_HighVsLowGainTS3.insert( std::pair<int,TH2F*>(ib*10+iski, htmp2) );   

      htmp2=gaindir.make<TH2F>("LowGainVsTOTTS3","LowGainVsTOTTS3",150, 0.,1500,175,-500.,3000);
      m_h_LowGainVsTOTTS3.insert( std::pair<int,TH2F*>(ib*10+iski, htmp2) );   
      
    }
  }
  m_h_nhitUnderSat=fs->make<TH1F>("Nhit_With_Negative_Saturation","Nhit_With_Negative_Saturation",1200,0,1200);
  m_h_nhitUnderSat_VS_SkirocID=fs->make<TH2F>("Nhit_With_Negative_Saturation_VS_SkirocID","Nhit_With_Negative_Saturation_VS_SkirocID",40,0,40,32,0,32);
  m_h_EEQuality=fs->make<TH1F>("EEQuality","EEQuality",3,0,3);
  m_h_FHQuality=fs->make<TH1F>("FHQuality","FHQuality",9,0,9);
  m_h_CenterQuality=fs->make<TH1F>("CenterQuality","CenterQuality",7,0,7);
  m_h_FullQuality=fs->make<TH1F>("FullQuality","FullQuality",10,0,10);
  m_tree=fs->make<TTree>("tree","Tree with common mode noise (one entry per event per skiroc)");
  m_tree->Branch("eventID",&m_evtID);
  m_tree->Branch("layerID",&m_layerID);
  m_tree->Branch("skirocID",&m_skirocID);
  m_tree->Branch("cmHighGain",&m_cmHigh,"cmHighGain[11]/F");
  m_tree->Branch("cmLowGain",&m_cmLow,"cmLowGain[11]/F");
  m_tree->Branch("time",&m_time,"time[11]/I");
}

void RawHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  if( !hits->size() )
    return;

  std::map<int,TH2Poly*>  polyMap;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    int Board_IU = 0;
    int Board_IV = 0; 
    int Board_Layer = 0;
    for(int ib = 0; ib<m_NHexaBoards; ib++) {
      Board_IU = 0;
      Board_IV = 0; 
      if( (ib == 6) || (ib == 9) ){
  Board_IU = 0;
  Board_IV = -1;
      }
      if( (ib == 5) || (ib == 8) ){
        Board_IU = 1;
        Board_IV = -1;
      }
      if(ib <= 3) Board_Layer = ib + 1;
      else if( (ib == 4) || (ib == 5) || (ib == 6) ) Board_Layer = 5;
      else if( (ib == 7) || (ib == 8) || (ib == 9) ) Board_Layer = 6;
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
  TH2Poly *h=dir.make<TH2Poly>();
  os.str("");
  os<<"HexaBoard"<<ib<<"_TimeSample"<<it;
  h->SetName(os.str().c_str());
  h->SetTitle(os.str().c_str());
  InitTH2Poly(*h, Board_Layer, Board_IU, Board_IV);
  polyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      }
    }
  }
  
  m_evtID+=1;
  CommonMode cm(essource_.emap_, true, true, -1.);
  cm.Evaluate( hits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();
  for( std::map<int,commonModeNoise>::iterator it=cmMap.begin(); it!=cmMap.end(); ++it ){
    m_skirocID=(it->first)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    m_layerID=(it->first)/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int key=m_layerID*1000+m_skirocID*100;

    for( uint16_t ts=0; ts<NUMBER_OF_TIME_SAMPLES; ts++ ){
      m_cmHigh[ts]=it->second.fullHG[ts];
      m_cmLow[ts]=it->second.fullLG[ts];
      m_time[ts]=ts;
      m_h_cmHigh[key]->Fill( m_cmHigh[ts] );
      m_h_cmLow[key]->Fill( m_cmLow[ts] );
      if( m_cmHigh[ts]>-160 && m_cmHigh[0]!=0 )
	m_h_cmHighTime[m_layerID*10+m_skirocID]->Fill( ts,(m_cmHigh[ts]-m_cmHigh[0])/m_cmHigh[0] );
      if( m_cmLow[ts]>-160 && m_cmLow[0]!=0 )
      m_h_cmLowTime[m_layerID*10+m_skirocID]->Fill( ts,(m_cmLow[ts]-m_cmLow[0])/m_cmLow[0] );
      key+=1;
    }
    m_tree->Fill();
  }

  int nhitUnderSat[m_NHexaBoards*HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA];
  for( int i=0; i<m_NHexaBoards*HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; i++ ) nhitUnderSat[i]=0;
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iboard=hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ichan=hit.channel();
    int iski=hit.skiroc();
    
    m_h_HighVsLowGainTS3[10*iboard+(iski%4)]->Fill(hit.lowGainADC(3), hit.highGainADC(3));
    m_h_LowGainVsTOTTS3[10*iboard+(iski%4)]->Fill(hit.totSlow(), hit.lowGainADC(3));

    std::pair<int,HGCalTBDetId> p( iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan,hit.detid() );
    setOfConnectedDetId.insert(p);
    channelInfo* cif=m_channelMap[iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan];
    if( hit.isUnderSaturationForHighGain() ) nhitUnderSat[hit.skiroc()]++;
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float highGain,lowGain;
      if( m_subtractCommonMode ){
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
      polyMap[ 100*iboard+it ]->Fill(CellCentreXY.first , CellCentreXY.second, highGain);
    }
  }
  int ntot=0;
  int nhitUnderZeroPerLayer[m_NHexaBoards];
  for( int i=0; i<m_NHexaBoards; i++ ) nhitUnderZeroPerLayer[i]=0;
  for( int i=0; i<m_NHexaBoards*HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; i++ ) {
    int iboard=i/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int iskiroc=i%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    ntot+=nhitUnderSat[i];
    m_h_nhitUnderSat_VS_SkirocID->Fill(i,nhitUnderSat[i]);
    m_h_nhitPerSkirocUnderSat[10*iboard+iskiroc]->Fill(nhitUnderSat[i]);
    nhitUnderZeroPerLayer[iboard]+=nhitUnderSat[i];
  }
  m_h_nhitUnderSat->Fill(ntot);
  int eeq(0),fhq(0),centerq(0),fullq(0);
  for( int i=0; i<m_NHexaBoards; i++ ){
    if( nhitUnderZeroPerLayer[i]<32 )
      continue;
    fullq++;
    if( i==0 || i==1 ){
      eeq++;
      centerq++;
      continue;
    }
    else
      fhq++;
    if( i<5||i==7 )
      centerq++;
  }
  m_h_EEQuality->Fill(eeq);
  m_h_FHQuality->Fill(fhq);
  m_h_CenterQuality->Fill(centerq);
  m_h_FullQuality->Fill(fullq);
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
  int Board_IU = 0;
  int Board_IV = 0; 
  int Board_Layer = 0;
  for(int ib = 0; ib<m_NHexaBoards; ib++) {
    Board_IU = 0;
    Board_IV = 0; 
    if( (ib == 6) || (ib == 9) ){
  Board_IU = 0;
  Board_IV = -1;
    }
    if( (ib == 5) || (ib == 8) ){
        Board_IU = 1;
        Board_IV = -1;
    }
    if(ib <= 3) Board_Layer = ib + 1;
    else if( (ib == 4) || (ib == 5) || (ib == 6) ) Board_Layer = 5;
    else if( (ib == 7) || (ib == 8) || (ib == 9) ) Board_Layer = 6;
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
      InitTH2Poly(*h, Board_Layer, Board_IU, Board_IV);
      hgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgpdir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, Board_Layer, Board_IU, Board_IV);
      lgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=hgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, Board_Layer, Board_IU, Board_IV);
      hgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, Board_Layer, Board_IU, Board_IV);
      lgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
    }
  }
  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = CellCentreXY.first;
    double iuy = CellCentreXY.second;
    channelInfo* cif=m_channelMap[iboard*1000+iski*100+ichan];
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float hgMean=cif->meanHGMap[it]/cif->counterHGMap[it];
      float lgMean=cif->meanLGMap[it]/cif->counterLGMap[it];
      float hgRMS=std::sqrt(cif->rmsHGMap[it]/cif->counterHGMap[it]-cif->meanHGMap[it]/cif->counterHGMap[it]*cif->meanHGMap[it]/cif->counterHGMap[it]);
      float lgRMS=std::sqrt(cif->rmsLGMap[it]/cif->counterLGMap[it]-cif->meanLGMap[it]/cif->counterLGMap[it]*cif->meanLGMap[it]/cif->counterLGMap[it]);
      hgMeanMap[ 100*iboard+it ]->Fill(iux , iuy, hgMean );
      lgMeanMap[ 100*iboard+it ]->Fill(iux , iuy, lgMean );
      hgRMSMap[ 100*iboard+it ]->Fill(iux , iuy, hgRMS );
      lgRMSMap[ 100*iboard+it ]->Fill(iux , iuy, lgRMS );
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
