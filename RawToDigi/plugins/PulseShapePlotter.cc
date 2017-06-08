#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include <fstream>
#include <sstream>
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

#include <iomanip>
#include <set>

const static size_t N_HEXABOARDS = 1;
const static size_t N_SKIROC_PER_HEXA = 4;
const static size_t N_CHANNELS_PER_SKIROC = 64;

class PulseShapePlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit PulseShapePlotter(const edm::ParameterSet&);
  ~PulseShapePlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  std::string m_electronicMap;
  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  double m_commonModeThreshold; //number of sigmas from ped mean

  int m_evtID;
  edm::EDGetTokenT<HGCalTBRawHitCollection> m_HGCalTBRawHitCollection;

  std::string m_pedestalHigh_filename;
  std::string m_pedestalLow_filename;

  struct pedestalChannel{
    HGCalTBDetId id;
    float pedHGMean[NUMBER_OF_TIME_SAMPLES];
    float pedLGMean[NUMBER_OF_TIME_SAMPLES];
    float pedHGRMS[NUMBER_OF_TIME_SAMPLES];
    float pedLGRMS[NUMBER_OF_TIME_SAMPLES];
  };
  std::map<int,pedestalChannel> m_pedMap; //key=10000*hexa+100*chip+chan

  struct commonModeNoise{
    commonModeNoise():fullHG(0),halfHG(0),mouseBiteHG(0),outerHG(0),fullLG(0),halfLG(0),mouseBiteLG(0),outerLG(0),fullCounter(0),halfCounter(0),mouseBiteCounter(0),outerCounter(0){;}
    float fullHG,halfHG,mouseBiteHG,outerHG;
    float fullLG,halfLG,mouseBiteLG,outerLG;
    int fullCounter,halfCounter,mouseBiteCounter,outerCounter;
  };

};

PulseShapePlotter::PulseShapePlotter(const edm::ParameterSet& iConfig) :
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneLayers_May2017.txt")),
  m_commonModeThreshold(iConfig.getUntrackedParameter<double>("CommonModeThreshold",3)),
  m_pedestalHigh_filename(iConfig.getParameter<std::string>("HighGainPedestalFileName")),
  m_pedestalLow_filename(iConfig.getParameter<std::string>("LowGainPedestalFileName"))
{
  m_HGCalTBRawHitCollection = consumes<HGCalTBRawHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));
  m_evtID=0;
  std::cout << iConfig.dump() << std::endl;
}


PulseShapePlotter::~PulseShapePlotter()
{

}

void PulseShapePlotter::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  FILE* file;
  char buffer[300];
  file = fopen (m_pedestalHigh_filename.c_str() , "r");
  if (file == NULL){ perror ("Error opening pedestal high gain"); exit(1); }
  else{
    while ( ! feof (file) ){
      if ( fgets (buffer , 300 , file) == NULL ) break;
      pedestalChannel ped;
      const char* index = buffer;
      int hexaboard,skiroc,channel,ptr,nval;
      nval=sscanf( index, "%d %d %d %n",&hexaboard,&skiroc,&channel,&ptr );
      if( nval==3 ){
	HGCalTBElectronicsId eid;
	switch( skiroc ){
	case 0 : eid=HGCalTBElectronicsId( 1, channel);break;
	case 1 : eid=HGCalTBElectronicsId( 4, channel);break;
	case 2 : eid=HGCalTBElectronicsId( 3, channel);break;
	case 3 : eid=HGCalTBElectronicsId( 2, channel);break;
	}
	if (!essource_.emap_.existsEId(eid.rawId()))
	  ped.id = HGCalTBDetId(-1);
	else
	  ped.id = essource_.emap_.eid2detId(eid);
	index+=ptr;
      }else continue;
      for( unsigned int ii=0; ii<NUMBER_OF_TIME_SAMPLES-0; ii++ ){
	float mean,rms;
	nval = sscanf( index, "%f %f %n",&mean,&rms,&ptr );
	if( nval==2 ){
	  ped.pedHGMean[ii]=mean;
	  ped.pedHGRMS[ii]=rms;
	  index+=ptr;
	}else continue;
      }
      m_pedMap.insert( std::pair<int,pedestalChannel>(10000*hexaboard+100*skiroc+channel,ped) );
    }
    fclose (file);
  }
  
  file = fopen (m_pedestalLow_filename.c_str() , "r");
  if (file == NULL){ perror ("Error opening pedestal low gain"); exit(1); }
  else{
    while ( ! feof (file) ){
      if ( fgets (buffer , 300 , file) == NULL ) break;
      pedestalChannel ped;
      const char* index = buffer;
      int hexaboard,skiroc,channel,ptr,nval,key;
      nval=sscanf( index, "%d %d %d %n",&hexaboard,&skiroc,&channel,&ptr );
      if( nval==3 ){
	key=10000*hexaboard+100*skiroc+channel;
	index+=ptr;
      }else continue;
      for( unsigned int ii=0; ii<NUMBER_OF_TIME_SAMPLES-0; ii++ ){
	float mean,rms;
	nval = sscanf( index, "%f %f %n",&mean,&rms,&ptr );
	if( nval==2 ){
	  m_pedMap[key].pedLGMean[ii]=mean;
	  m_pedMap[key].pedLGRMS[ii]=rms;
	  index+=ptr;
	}else continue;
      }
    }
    fclose (file);
  }
}

void PulseShapePlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  std::map<int,TH1F*>  hMapHG,hMapLG;
  std::ostringstream os( std::ostringstream::ate );
  os << "Event" << event.id().event();
  TFileDirectory dir = fs->mkdir( os.str().c_str() );
  for(size_t ib = 0; ib<N_HEXABOARDS; ib++) {
    for(size_t is = 0; is<N_SKIROC_PER_HEXA; is++) {
      for( size_t ic=0; ic<N_CHANNELS_PER_SKIROC; ic++ ){
	HGCalTBElectronicsId eid;
	switch( is ){
	case 0 : eid=HGCalTBElectronicsId( 1, ic);break;
	case 1 : eid=HGCalTBElectronicsId( 4, ic);break;
	case 2 : eid=HGCalTBElectronicsId( 3, ic);break;
	case 3 : eid=HGCalTBElectronicsId( 2, ic);break;
	}
	if (!essource_.emap_.existsEId(eid.rawId())) continue;
	os.str("");
	os<<"HexaBoard"<<ib<<"_Skiroc"<<is<<"_Channel"<<ic<<"_HG";
	TH1F *hHG=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES,0,NUMBER_OF_TIME_SAMPLES*25);
	hMapHG.insert( std::pair<int,TH1F*>(1000*ib+100*is+ic,hHG) );
	os.str("");
	os<<"HexaBoard"<<ib<<"_Skiroc"<<is<<"_Channel"<<ic<<"_LG";
	TH1F* hLG=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),NUMBER_OF_TIME_SAMPLES,0,NUMBER_OF_TIME_SAMPLES*25);
	hMapLG.insert( std::pair<int,TH1F*>(1000*ib+100*is+ic,hLG) );
      }
    }
  }
  
  edm::Handle<HGCalTBRawHitCollection> hits;
  event.getByToken(m_HGCalTBRawHitCollection, hits);
  
  commonModeNoise cm[NUMBER_OF_TIME_SAMPLES-0][4];
  for( auto hit : *hits ){
    int iboard=hit.skiroc()/N_SKIROC_PER_HEXA;
    int iski=hit.skiroc();
    int ichan=hit.channel();
    if( !essource_.emap_.existsDetId(hit.detid()) ) continue;
    for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
      if( fabs(hit.highGainADC(it)-m_pedMap[10000*iboard+100*iski+ichan].pedHGMean[it])>m_commonModeThreshold*m_pedMap[10000*iboard+100*iski+ichan].pedHGRMS[it] ) continue;
      float highGain = hit.highGainADC(it)-m_pedMap[10000*iboard+100*iski+ichan].pedHGMean[it];
      float lowGain = hit.lowGainADC(it)-m_pedMap[10000*iboard+100*iski+ichan].pedLGMean[it];
      switch ( hit.detid().cellType() ){
      case 0 : cm[it][iski].fullHG += highGain; cm[it][iski].fullLG += lowGain; cm[it][iski].fullCounter++; break;
      case 2 : cm[it][iski].halfHG += highGain; cm[it][iski].halfLG += lowGain; cm[it][iski].halfCounter++; break;
      case 3 : cm[it][iski].mouseBiteHG += highGain; cm[it][iski].mouseBiteLG += lowGain; cm[it][iski].mouseBiteCounter++; break;
      case 4 : cm[it][iski].outerHG += highGain; cm[it][iski].outerLG += lowGain; cm[it][iski].outerCounter++; break;
      }
    }
  }

  for( auto hit : *hits ){
    int iboard=hit.skiroc()/N_SKIROC_PER_HEXA;
    int iski=hit.skiroc();
    int ichan=hit.channel();
    if( essource_.emap_.existsDetId(hit.detid()) ){
      for( size_t it=0; it<NUMBER_OF_TIME_SAMPLES; it++ ){
	float highGain,lowGain;
	float subHG(0),subLG(0);
	switch ( hit.detid().cellType() ){
	case 0 : subHG=cm[it][iski].fullHG/cm[it][iski].fullCounter; subLG=cm[it][iski].fullLG/cm[it][iski].fullCounter; break;
	case 2 : subHG=cm[it][iski].halfHG/cm[it][iski].halfCounter; subLG=cm[it][iski].halfLG/cm[it][iski].halfCounter; break;
	case 3 : subHG=cm[it][iski].mouseBiteHG/cm[it][iski].mouseBiteCounter; subLG=cm[it][iski].mouseBiteLG/cm[it][iski].mouseBiteCounter; break;
	case 4 : subHG=cm[it][iski].outerHG/cm[it][iski].outerCounter; subLG=cm[it][iski].outerLG/cm[it][iski].outerCounter; break;
	}
	highGain=hit.highGainADC(it)-m_pedMap[10000*iboard+100*iski+ichan].pedHGMean[it]-subHG;
	lowGain=hit.lowGainADC(it)-m_pedMap[10000*iboard+100*iski+ichan].pedLGMean[it]-subLG;
	hMapHG[1000*iboard+100*iski+ichan]->Fill(25*it,highGain);
	hMapLG[1000*iboard+100*iski+ichan]->Fill(25*it,lowGain);
      }
    }
  }
  m_evtID++;
}

void PulseShapePlotter::endJob()
{
}

void PulseShapePlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PulseShapePlotter);
