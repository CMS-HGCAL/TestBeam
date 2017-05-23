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
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include <iomanip>
#include <set>

const static size_t N_HEXABOARDS = 1;
const static size_t N_TIME_SAMPLES = 13;
const static size_t N_SKIROC_PER_HEXA = 4;
const static size_t N_CHANNELS_PER_SKIROC = 64;

#define MAXVERTICES 6
static const double delta = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.

class RawDataPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit RawDataPlotter(const edm::ParameterSet&);
  ~RawDataPlotter();
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
  int m_sensorsize;
  bool m_eventPlotter;

  int m_evtID;
  std::map<int,TH1F*> m_h_adcHigh;
  std::map<int,TH1F*> m_h_adcLow;
  std::map<int,TH2F*> m_h_pulseHigh;
  std::map<int,TH2F*> m_h_pulseLow;

  edm::EDGetTokenT<HGCalTBSkiroc2CMSCollection> m_HGCalTBSkiroc2CMSCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;
  std::string m_pedHighGainFileName;
  std::string m_pedLowGainFileName;

};

RawDataPlotter::RawDataPlotter(const edm::ParameterSet& iConfig) :
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
  m_pedHighGainFileName(iConfig.getParameter<std::string>("HighGainPedestalFileName")),
  m_pedLowGainFileName(iConfig.getParameter<std::string>("LowGainPedestalFileName"))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_evtID=0;
  
  std::ostringstream os( std::ostringstream::ate );
  TH2F* htmp2;
  TH1F* htmp1;
  for(size_t ib = 0; ib<N_HEXABOARDS; ib++) {
    for( size_t iski=0; iski<N_SKIROC_PER_HEXA; iski++ ){
      os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
      TFileDirectory dir = fs->mkdir( os.str().c_str() );
      for( size_t ichan=0; ichan<N_CHANNELS_PER_SKIROC; ichan++ ){
	for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
	  os.str("");
	  os << "HighGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_Sample" << it ;
	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
	  m_h_adcHigh.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
	  os.str("");
	  os << "LowGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_Sample" << it ;
	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
	  m_h_adcLow.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
	}
	if( ichan%2==0 ){
	  os.str("");
	  os << "PulseHighGain_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	  htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),N_TIME_SAMPLES-2,0, (N_TIME_SAMPLES-2)*25,1000,0,4096);
	  m_h_pulseHigh.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp2) );
	  os.str("");
	  os << "PulseLowGain_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	  htmp2=dir.make<TH2F>(os.str().c_str(),os.str().c_str(),N_TIME_SAMPLES-2,0, (N_TIME_SAMPLES-2)*25,1000,0,4096);
	  m_h_pulseLow.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp2) );
	}
      }
    }
  }
  std::cout << iConfig.dump() << std::endl;
}


RawDataPlotter::~RawDataPlotter()
{

}

void RawDataPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);
  
  std::map<int,TH2Poly*>  polyMap;
  if( m_eventPlotter ){
    std::ostringstream os( std::ostringstream::ate );
    os << "Event" << event.id().event();
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    for(size_t ib = 0; ib<N_HEXABOARDS; ib++) {
      for( size_t it=0; it<N_TIME_SAMPLES-2; it++ ){
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
  
  for( size_t iski=0;iski<skirocs->size(); iski++ ){
    HGCalTBSkiroc2CMS skiroc=skirocs->at(iski);
    std::vector<int> rollpositions=skiroc.rollPositions();
    int iboard=iski/N_SKIROC_PER_HEXA;
    for( size_t ichan=0; ichan<N_CHANNELS_PER_SKIROC; ichan++ ){
      for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
	if( rollpositions[it]!=11 && rollpositions[it]!=12 ){ //rm on track samples
	  m_h_adcHigh[iboard*100000+iski*10000+ichan*100+it]->Fill(skiroc.ADCHigh(ichan,it));
	  m_h_adcLow[iboard*100000+iski*10000+ichan*100+it]->Fill(skiroc.ADCLow(ichan,it));
	  if( ichan%2==0 ){
	    m_h_pulseHigh[iboard*1000+iski*100+ichan]->Fill( rollpositions[it]*25,skiroc.ADCHigh(ichan,it) ); 
	    m_h_pulseLow[iboard*1000+iski*100+ichan]->Fill( rollpositions[it]*25,skiroc.ADCLow(ichan,it) );
	    
	    HGCalTBDetId detid=skiroc.detid( ichan );
	    if(!IsCellValid.iu_iv_valid( detid.layer(),detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize ) )  continue;
	    if(m_eventPlotter){
	      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
	      double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
	      double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
	      polyMap[ 100*iboard+rollpositions[it] ]->Fill(iux/2 , iuy, skiroc.ADCHigh(ichan,it) );
	    }
	    std::pair<int,HGCalTBDetId> p( iboard*1000+iski*100+ichan,detid );
	    setOfConnectedDetId.insert(p);
	  }
	}
      }
    }
  }
}

void RawDataPlotter::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
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

void RawDataPlotter::beginJob()
{
}

// const static size_t N_HEXABOARDS = 1;
// const static size_t N_TIME_SAMPLES = 13;
// const static size_t N_SKIROC_PER_HEXA = 4;
// const static size_t N_CHANNELS_PER_SKIROC = 64;
void RawDataPlotter::endJob()
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  TFileDirectory dir = fs->mkdir( "PedestalPlotter" );
  std::map<int,TH2Poly*>  pedPolyMap;
  std::map<int,TH2Poly*>  pedPolyMapNC;
  std::map<int,TH2Poly*>  noisePolyMap;
  std::map<int,TH2Poly*>  noisePolyMapNC;
  std::map<int,TH2Poly*>  chanMap;
  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h;
  for(size_t ib = 0; ib<N_HEXABOARDS; ib++) {
    for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"HighGain_HexaBoard"<<ib<<"_SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      pedPolyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"NC_HighGain_HexaBoard"<<ib<<"_SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      pedPolyMapNC.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"Noise_HighGain_HexaBoard"<<ib<<"_SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      noisePolyMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
      h=dir.make<TH2Poly>();
      os.str("");
      os<<"NC_Noise_HighGain_HexaBoard"<<ib<<"_SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      noisePolyMapNC.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
    }
    h=dir.make<TH2Poly>();
    os.str("");
    os<<"ChannelMapping";
    h->SetName(os.str().c_str());
    h->SetTitle(os.str().c_str());
    InitTH2Poly(*h, (int)ib, 0, 0);
    chanMap.insert( std::pair<int,TH2Poly*>(ib,h) );
  }

  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    int ichanNC=(*it).first%100+1;
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + delta) : (CellCentreXY.first - delta) ;
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + delta) : (CellCentreXY.second - delta);
    for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
      pedPolyMap[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcHigh[iboard*100000+iski*10000+ichan*100+it]->GetMean() );
      pedPolyMapNC[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcHigh[iboard*100000+iski*10000+ichanNC*100+it]->GetMean() );
      noisePolyMap[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcHigh[iboard*100000+iski*10000+ichan*100+it]->GetRMS() );
      noisePolyMapNC[ 100*iboard+it ]->Fill(iux/2 , iuy, m_h_adcHigh[iboard*100000+iski*10000+ichanNC*100+it]->GetRMS() );
    }
    chanMap[ iboard ]->Fill(iux/2 , iuy, iski*1000+ichan );
  }

  std::fstream pedestalHG;pedestalHG.open(m_pedHighGainFileName,std::ios::out);
  std::fstream pedestalLG;pedestalLG.open(m_pedLowGainFileName,std::ios::out);
  for(size_t ib = 0; ib<N_HEXABOARDS; ib++) {
    for( size_t iski=0; iski<N_SKIROC_PER_HEXA; iski++ ){
      for( size_t ichan=0; ichan<N_CHANNELS_PER_SKIROC; ichan++ ){
	pedestalHG << ib << " " << iski << " " << ichan << " " ;
	pedestalLG << ib << " " << iski << " " << ichan << " " ;
	for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
	  int key=ib*100000+iski*10000+ichan*100+it;
	  pedestalHG << m_h_adcHigh[key]->GetMean() << " " << m_h_adcHigh[key]->GetRMS() << " ";
	  pedestalLG << m_h_adcLow[key]->GetMean() << " " << m_h_adcLow[key]->GetRMS() << " ";
	}
	pedestalHG << "\n" ;
	pedestalLG << "\n" ;	
      }
    }
  }
  pedestalHG.close();
  pedestalLG.close();
}

void RawDataPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RawDataPlotter);
