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
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBLayer.h"
#include "HGCal/DataFormats/interface/HGCalTBModule.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Reco/interface/CommonMode.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
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
};

class chipInfo{
public:
  chipInfo()
  {
    for(int i=0; i<NUMBER_OF_TIME_SAMPLES; i++){
      hgsum.push_back(0);
      lgsum.push_back(0);
    }
  }
  void reset()
  {
    for(int i=0; i<NUMBER_OF_TIME_SAMPLES; i++){
      hgsum[i]=0;
      lgsum[i]=0;
    }
  }
  int chip,board,module;
  commonModeNoise cmn;
  std::vector<float> hgsum;
  std::vector<float> lgsum;
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
  void InitTH2Poly(TH2Poly& poly, int det, int layerID);//detector=0 for ee, =1 for fh

  std::string m_electronicMap;
  std::string m_detectorLayoutFile;
  HGCalElectronicsMap m_emap;
  HGCalTBDetectorLayout m_layout;

  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  int m_sensorsize;
  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;

  std::map<int,channelInfo*> m_channelMap;
  std::map<int,chipInfo> m_chipMap;

  TTree* m_tree;
  int m_evtID;
  int m_board;
  int m_chip;
  int m_module;
  float m_cmhg[NUMBER_OF_TIME_SAMPLES];
  float m_cmlg[NUMBER_OF_TIME_SAMPLES];
  float m_hgSum[NUMBER_OF_TIME_SAMPLES];
  float m_lgSum[NUMBER_OF_TIME_SAMPLES];
};

RawHitPlotter::RawHitPlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_detectorLayoutFile(iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128))
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
  if (!io.load(fip.fullPath(), m_emap)) {
    throw cms::Exception("Unable to load electronics map");
  };
  fip=edm::FileInPath(m_detectorLayoutFile);
  if (!io.load(fip.fullPath(), m_layout)) {
    throw cms::Exception("Unable to load detector layout file");
  };
  for( auto layer : m_layout.layers() )
    layer.print();
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  std::ostringstream os( std::ostringstream::ate );
  for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
    for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){ 
      int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*ib+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iski)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
      for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
	HGCalTBElectronicsId eid(skiId,ichan);      
	if( !m_emap.existsEId(eid) ) continue;
	int key=ib*1000+iski*100+ichan;
	channelInfo *cif=new channelInfo();
	cif->key=key;
	cif->init();
	m_channelMap.insert( std::pair<int,channelInfo*>(key,cif) );
      }
    }
  }
  m_tree=fs->make<TTree>("tree","Tree with common mode noise (one entry per event per skiroc)");
  m_tree->Branch("eventID",&m_evtID);
  m_tree->Branch("chip",&m_chip);
  m_tree->Branch("board",&m_board);
  m_tree->Branch("module",&m_module);
  m_tree->Branch("cmhg",&m_cmhg,"cmhg[11]/F");
  m_tree->Branch("cmlg",&m_cmlg,"cmlg[11]/F");
  m_tree->Branch("hgSum",&m_hgSum,"hgSum[11]/F");
  m_tree->Branch("lgSum",&m_lgSum,"lgSum[11]/F");
}

void RawHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  if( !hits->size() )
    return;

  CommonMode cm(m_emap); //default is common mode per chip using the median
  cm.Evaluate( hits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();

  for( std::map<int,chipInfo>::iterator it=m_chipMap.begin(); it!=m_chipMap.end(); ++it )
    it->second.reset();
  
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( m_emap.detId2eid(hit.detid().rawId()) );
    if( !m_emap.existsEId(eid) ) continue;
    int board=hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int channel=hit.channel();
    int chip=hit.skiroc(); //from 0 to NHexaboard*4-1
    HGCalTBLayer alayer = m_layout.at( hit.detid().layer()-1 );
    HGCalTBModule amodule = alayer.at( hit.detid().sensorIU(), hit.detid().sensorIV() );
    int module=amodule.moduleID();
    std::pair<int,HGCalTBDetId> p( board*1000+(chip%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+channel,hit.detid() );
    setOfConnectedDetId.insert(p);
    int chipID=chip%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;//from 0 to 3
    if( m_chipMap.find(module*100+chipID)==m_chipMap.end() ) {
      chipInfo c;
      c.chip=chipID;
      c.board=board;
      c.module=module;
      m_chipMap.insert( std::pair<int,chipInfo>(module*100+c.chip,c) );
    }
    m_chipMap[module*100+chipID].cmn=cmMap[chip];
    channelInfo* cif=m_channelMap[board*1000+(chip%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+channel];
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float highGain,lowGain;
      float subHG(0),subLG(0);
      switch ( hit.detid().cellType() ){
      case 0 : subHG=cmMap[chip].fullHG[it]; subLG=cmMap[chip].fullLG[it]; break;
      case 2 : subHG=cmMap[chip].halfHG[it]; subLG=cmMap[chip].halfLG[it]; break;
      case 3 : subHG=cmMap[chip].mouseBiteHG[it]; subLG=cmMap[chip].mouseBiteLG[it]; break;
      case 4 : subHG=cmMap[chip].outerHG[it]; subLG=cmMap[chip].outerLG[it]; break;
      case 5 : subHG=cmMap[chip].mergedHG[it]; subLG=cmMap[chip].mergedLG[it]; break;
      }
      highGain=hit.highGainADC(it)-subHG;
      lowGain=hit.lowGainADC(it)-subLG;

      cif->meanHGMap[it]+=highGain;
      cif->rmsHGMap[it]+=highGain*highGain;
      cif->counterHGMap[it]+=1;
      cif->meanLGMap[it]+=lowGain;
      cif->rmsLGMap[it]+=lowGain*lowGain;
      cif->counterLGMap[it]+=1;
      m_chipMap[module*100+chipID].hgsum[it]+=highGain;
      m_chipMap[module*100+chipID].lgsum[it]+=lowGain;
    }
  }
  for( std::map<int,chipInfo>::iterator it=m_chipMap.begin(); it!=m_chipMap.end(); ++it ){
    if( m_evtID==0 )
      std::cout << "End of 1st event:\n" << it->first << " " << it->second.module << " " << it->second.board << " " << it->second.chip << std::endl;
    m_chip=it->second.chip;
    m_board=it->second.board;
    m_module=it->second.module;
    for(int time=0; time<NUMBER_OF_TIME_SAMPLES; time++){
      m_cmhg[time]=it->second.cmn.fullHG[time];
      m_cmlg[time]=it->second.cmn.fullLG[time];
      m_hgSum[time]=it->second.hgsum[time];
      m_lgSum[time]=it->second.lgsum[time];
    }
    m_tree->Fill();
  }
  m_evtID+=1;
}

void RawHitPlotter::InitTH2Poly(TH2Poly& poly, int det, int layerID)
{
  double HexX[MAXVERTICES] = {0.};
  double HexY[MAXVERTICES] = {0.};
  if( det==0 ){
    for(int iv = -7; iv < 8; iv++) {
      for(int iu = -7; iu < 8; iu++) {
	if(!IsCellValid.iu_iv_valid(layerID, 0, 0, iu, iv, m_sensorsize)) 
	  continue;
	CellXY = TheCell.GetCellCoordinatesForPlots(layerID, 0, 0, iu, iv, m_sensorsize);
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
  else if( det==1 ){
    for(int sensorIV = -1; sensorIV <= 1; sensorIV++){
      for(int sensorIU = -1; sensorIU <= 1; sensorIU++){
	for(int iv = -7; iv < 8; iv++) {
	  for(int iu = -7; iu < 8; iu++) {
	    if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, m_sensorsize)) 
	      continue;
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
    }// loop over Sensor_Iu ends here
  }// loop over Sensor_Iv ends here
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
  for(int il = 0; il<m_layout.nlayers(); il++) {
    os.str("");
    os<<"Layer"<<il;
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    TFileDirectory hgpdir = dir.mkdir( "HighGainPedestal" );
    TFileDirectory lgpdir = dir.mkdir( "LowGainPedestal" );
    TFileDirectory hgndir = dir.mkdir( "HighGainNoise" );
    TFileDirectory lgndir = dir.mkdir( "LowGainNoise" );
    int subdetId = m_layout.at(il).subdet();
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      h=hgpdir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, subdetId, il);
      hgMeanMap.insert( std::pair<int,TH2Poly*>(100*il+it,h) );

      h=lgpdir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, subdetId, il);
      lgMeanMap.insert( std::pair<int,TH2Poly*>(100*il+it,h) );

      h=hgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, subdetId, il);
      hgRMSMap.insert( std::pair<int,TH2Poly*>(100*il+it,h) );

      h=lgndir.make<TH2Poly>();
      os.str("");
      os<<"TS"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, subdetId, il);
      lgRMSMap.insert( std::pair<int,TH2Poly*>(100*il+it,h) );
    }
  }
  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    HGCalTBDetId detid=(*it).second;
    int ilayer=detid.layer()-1;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = CellCentreXY.first;
    double iuy = CellCentreXY.second;
    channelInfo* cif=m_channelMap[iboard*1000+iski*100+ichan];
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      float hgMean=cif->meanHGMap[it]/cif->counterHGMap[it];
      float lgMean=cif->meanLGMap[it]/cif->counterLGMap[it];
      float hgRMS=std::sqrt(cif->rmsHGMap[it]/cif->counterHGMap[it]-cif->meanHGMap[it]/cif->counterHGMap[it]*cif->meanHGMap[it]/cif->counterHGMap[it]);
      float lgRMS=std::sqrt(cif->rmsLGMap[it]/cif->counterLGMap[it]-cif->meanLGMap[it]/cif->counterLGMap[it]*cif->meanLGMap[it]/cif->counterLGMap[it]);
      hgMeanMap[ 100*ilayer+it ]->Fill(iux , iuy, hgMean );
      lgMeanMap[ 100*ilayer+it ]->Fill(iux , iuy, lgMean );
      hgRMSMap[ 100*ilayer+it ]->Fill(iux , iuy, hgRMS );
      lgRMSMap[ 100*ilayer+it ]->Fill(iux , iuy, lgRMS );
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
