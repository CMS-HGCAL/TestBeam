#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
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

using namespace std;

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

class RawHit_EventDisplay : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit RawHit_EventDisplay(const edm::ParameterSet&);
  ~RawHit_EventDisplay();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTH2Poly(TH2Poly& poly, int layerID);

  std::string m_electronicMap;

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

  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;
  
};

RawHit_EventDisplay::RawHit_EventDisplay(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_July_6Layers.txt")),
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_subtractCommonMode(iConfig.getUntrackedParameter<bool>("SubtractCommonMode",false))
  {
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));
  
  m_evtID=0;
  
  std::cout << iConfig.dump() << std::endl;
}


RawHit_EventDisplay::~RawHit_EventDisplay()
{

}

void RawHit_EventDisplay::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  std::ostringstream os( std::ostringstream::ate );
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
	m_channelMap.insert( std::pair<int,channelInfo*>(key,cif) );
      }
    }
  }
}

void RawHit_EventDisplay::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  if( !hits->size() )
    return;

  if( (event.id().event() % 20) == 0 ){ // Plot event displays every 20 events(for now)

  std::map<int,TH2Poly*>  polyMap;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_LAYERS; ib++) {
	size_t it = 3;
//      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	TH2Poly *h=dir.make<TH2Poly>();
	os.str("");
	os<<"Layer"<<(ib+1)<<"_TimeSample"<<it;
	h->SetName(os.str().c_str());
	h->SetTitle(os.str().c_str());
	InitTH2Poly(*h, (int) (ib + 1));
	polyMap.insert( std::pair<int,TH2Poly*>(100*(ib + 1) + it, h) );
//      }// Loop over time samples

    }
  }
  
  CommonMode cm(essource_.emap_); //default is common mode per chip using the median
  cm.Evaluate( hits );
  std::map<int,commonModeNoise> cmMap=cm.CommonModeNoiseMap();
  
 
  for( auto hit : *hits ){
    HGCalTBElectronicsId eid( essource_.emap_.detId2eid(hit.detid().rawId()) );
    if( !essource_.emap_.existsEId(eid) ) continue;
    int iboard=hit.skiroc()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    int ichan=hit.channel();
    int iski=hit.skiroc();
//    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      size_t it = 3;    
      float highGain;
//      float lowGain;
      if( m_subtractCommonMode ){
  	iski = eid.iskiroc();
  	float subHG(0);
//	float subLG(0);
  	switch ( hit.detid().cellType() ){
  	case 0 :{
		 subHG=cmMap[iski].fullHG[it];
//		 subLG=cmMap[iski].fullLG[it];
		 break;
		}
  	case 2 :{
		 subHG=cmMap[iski].halfHG[it];
//		 subLG=cmMap[iski].halfLG[it];
		 break;
		}
  	case 3 :{
		 subHG=cmMap[iski].mouseBiteHG[it];
//		 subLG=cmMap[iski].mouseBiteLG[it];
		 break;
		}
  	case 4 :{
		 subHG=cmMap[iski].outerHG[it];
//		 subLG=cmMap[iski].outerLG[it];
		 break;
		}
  	}
  	highGain=hit.highGainADC(it)-subHG;
//  	lowGain=hit.lowGainADC(it)-subLG;
      }
      else{
  	highGain=hit.highGainADC(it);
//  	lowGain=hit.lowGainADC(it);
      }

      if(!m_eventPlotter||!IsCellValid.iu_iv_valid(hit.detid().layer(), hit.detid().sensorIU(), hit.detid().sensorIV(), hit.detid().iu(), hit.detid().iv(), m_sensorsize))  continue;
      CellCentreXY=TheCell.GetCellCentreCoordinatesForPlots(hit.detid().layer(), hit.detid().sensorIU(), hit.detid().sensorIV(), hit.detid().iu(), hit.detid().iv(), m_sensorsize);
      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
      polyMap[ 100*( hit.detid().layer() ) + it ]->Fill(iux , iuy, highGain);
//    } // Loop over time samples
    
  } // Loop over Hits

 }//IF CONDITION TO CONSIDER EVERY 20th EVENT

}// analyze method

void RawHit_EventDisplay::InitTH2Poly(TH2Poly& poly, int layerID){
double HexX[MAXVERTICES] = {0.};
double HexY[MAXVERTICES] = {0.};
	for(int SensorIV = -1; SensorIV <= 1; SensorIV++){
		for(int SensorIU = -1; SensorIU <= 1; SensorIU++){
			for(int iv = -7; iv < 8; iv++) {
				for(int iu = -7; iu < 8; iu++) {
					if(!IsCellValid.iu_iv_valid(layerID, SensorIU, SensorIV, iu, iv, m_sensorsize)) continue;
					CellXY = TheCell.GetCellCoordinatesForPlots(layerID, SensorIU, SensorIV, iu, iv, m_sensorsize);
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

		}// loop over Sensor_Iu ends here
	}// loop over Sensor_Iv ends here
}


void RawHit_EventDisplay::endJob()
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
}

void RawHit_EventDisplay::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RawHit_EventDisplay);
