#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TProfile.h"
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

const static size_t N_HEXABOARDS = 1;
const static size_t N_TIME_SAMPLES = 13;
const static size_t N_SKIROC_PER_HEXA = 4;
const static size_t N_CHANNELS_PER_SKIROC = 64;

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

  std::string m_electronicMap;

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  int m_sensorsize;

  TTree* m_tree;
  int m_evtID;
  std::map<int,float> m_adcHigh;
  std::map<int,float> m_adcLow;
  std::map<int,TH2F*> m_h_pulseHigh;
  std::map<int,TH2F*> m_h_pulseLow;

  edm::EDGetTokenT<HGCalTBSkiroc2CMSCollection> m_HGCalTBSkiroc2CMSCollection;
};

RawDataPlotter::RawDataPlotter(const edm::ParameterSet& iConfig) :
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize"))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  m_tree = fs->make<TTree>("tree", "HGCAL TB variables tree");

  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_evtID=0;
  
  m_tree->Branch( "evtID" , &m_evtID); 
  std::ostringstream os( std::ostringstream::ate );
  TH2F* htmp;
  for(size_t ib = 0; ib<N_HEXABOARDS; ib++) {
    for( size_t iski=0; iski<N_SKIROC_PER_HEXA; iski++ ){
      for( size_t ichan=0; ichan<N_CHANNELS_PER_SKIROC; ichan++ ){
	for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
	  m_adcHigh.insert( std::pair<int,float>(ib*100000+iski*10000+ichan*100+it, 0.) );
	  os.str("");
	  os << "HighGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_Sample" << it ;
	  m_tree->Branch( os.str().c_str(),&m_adcHigh[ib*100000+iski*10000+ichan*100+it]); 
	  m_adcLow.insert( std::pair<int,float>(ib*100000+iski*10000+ichan*100+it, 0.) );
	  os.str("");
	  os << "LowGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_Sample" << it ;
	  m_tree->Branch( os.str().c_str(),&m_adcLow[ib*100000+iski*10000+ichan*100+it]); 
	}
	if( ichan%2==0 ){
	  os.str("");
	  os << "PulseHighGain_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	  htmp=fs->make<TH2F>(os.str().c_str(),os.str().c_str(),N_TIME_SAMPLES-2,0, (N_TIME_SAMPLES-2)*25,1000,0,4096);
	  m_h_pulseHigh.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp) );
	  os.str("");
	  os << "PulseLowGain_Hexa" << ib << "_Chip" << iski << "_Channel" << ichan;
	  htmp=fs->make<TH2F>(os.str().c_str(),os.str().c_str(),N_TIME_SAMPLES-2,0, (N_TIME_SAMPLES-2)*25,1000,0,4096);
	  m_h_pulseLow.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp) );
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
  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);
  
  for( size_t iski=0;iski<skirocs->size(); iski++ ){
    HGCalTBSkiroc2CMS skiroc=skirocs->at(iski);
    std::vector<int> rollpositions=skiroc.rollPositions();
    int iboard=iski/N_SKIROC_PER_HEXA;
    for( size_t ichan=0; ichan<N_CHANNELS_PER_SKIROC; ichan++ ){
      for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
	if( rollpositions[it]!=11 && rollpositions[it]!=12 ){ //rm on track samples
	  m_adcHigh[iboard*100000+iski*10000+ichan*100+it]=skiroc.ADCHigh(ichan,it);
	  m_adcLow[iboard*100000+iski*10000+ichan*100+it]=skiroc.ADCLow(ichan,it);
	}
	else{
	  m_adcHigh[iboard*100000+iski*10000+ichan*100+it]=-1000;
	  m_adcLow[iboard*100000+iski*10000+ichan*100+it]=-1000;
	}
      }
      if( ichan%2==0 ){
	for( size_t it=0; it<N_TIME_SAMPLES; it++ ){
	  if( rollpositions[it]==11 || rollpositions[it]==12 ) continue;
	  else{
	    m_h_pulseHigh[iboard*1000+iski*100+ichan]->Fill( rollpositions[it]*25,skiroc.ADCHigh(ichan,it) ); 
	    m_h_pulseLow[iboard*1000+iski*100+ichan]->Fill( rollpositions[it]*25,skiroc.ADCLow(ichan,it) );
	  }
	}
      }
    }
  }
  m_tree->Fill();
}

void RawDataPlotter::beginJob()
{
}

void RawDataPlotter::endJob()
{
}

void RawDataPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RawDataPlotter);
