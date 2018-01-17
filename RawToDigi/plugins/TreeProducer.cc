#include <iostream>
#include "TTree.h"
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
#include <iomanip>
#include <set>

class TreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit TreeProducer(const edm::ParameterSet&);
  ~TreeProducer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  struct {
    HGCalElectronicsMap emap_;
  } essource_;

  bool m_subtractPedestal;
  std::string m_pedestalHigh_filename;
  std::string m_pedestalLow_filename;
  std::string m_electronicMap;
  int m_event;

  TTree* m_tree;
  int m_chip;
  int m_roll;
  int m_dacinj;
  int m_timesamp[13];
  int m_hg[13][64];
  int m_lg[13][64];
  int m_tot_fast[64];
  int m_tot_slow[64];
  int m_toa_rise[64];
  int m_toa_fall[64];

  edm::EDGetTokenT<HGCalTBSkiroc2CMSCollection> m_HGCalTBSkiroc2CMSCollection;
};

TreeProducer::TreeProducer(const edm::ParameterSet& iConfig) :
  m_subtractPedestal( iConfig.getUntrackedParameter<bool>("SubtractPedestal",false) ),
  m_pedestalHigh_filename( iConfig.getUntrackedParameter<std::string>("HighGainPedestalFileName",std::string("pedestalHG.txt")) ),
  m_pedestalLow_filename( iConfig.getUntrackedParameter<std::string>("LowGainPedestalFileName",std::string("pedestalLG.txt")) ),
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_OneModule.txt"))
{
  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_event=0;

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  std::cout << iConfig.dump() << std::endl;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  m_tree = fs->make<TTree>("sk2cms", "sk2cms tree");
  m_tree->Branch("event", &m_event);
  m_tree->Branch("chip", &m_chip);
  m_tree->Branch("roll", &m_roll);
  m_tree->Branch("dacinj", &m_dacinj);
  m_tree->Branch("timesamp", &m_timesamp,"timesamp[13]/I");
  m_tree->Branch("hg", &m_hg,"hg[13][64]/I");
  m_tree->Branch("lg", &m_lg,"lg[13][64]/I");
  m_tree->Branch("tot_fast", &m_tot_fast,"tot_fast[64]/I");
  m_tree->Branch("tot_slow", &m_tot_slow,"tot_slow[64]/I");
  m_tree->Branch("toa_rise", &m_toa_rise,"toa_rise[64]/I");
  m_tree->Branch("toa_fall", &m_toa_fall,"toa_fall[64]/I");

}


TreeProducer::~TreeProducer()
{

}

void TreeProducer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);

  if( !skirocs->size() ) return;
  
  for( size_t iski=0;iski<skirocs->size(); iski++ ){
    HGCalTBSkiroc2CMS skiroc=skirocs->at(iski);
    std::vector<int> rollpositions=skiroc.rollPositions();
    m_chip=iski;
    m_roll=skiroc.rollMask();
    m_dacinj=skiroc.dacInjection();
    std::vector<int> ts=skiroc.rollPositions();
    int its=0;
    for( std::vector<int>::iterator it=ts.begin(); it!=ts.end(); ++it ){
      m_timesamp[its]=(*it);
      its++;
    }
    for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
      //      if( ichan==10 )
      //	std::cout << std::dec << iski << " " << skiroc.TOTSlow(ichan) << std::endl;
      for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
	m_hg[it][ichan]=skiroc.ADCHigh(ichan,it);
	m_lg[it][ichan]=skiroc.ADCLow(ichan,it);
      }
      m_tot_fast[ichan]=skiroc.TOTFast(ichan);
      m_tot_slow[ichan]=skiroc.TOTSlow(ichan);
      m_toa_rise[ichan]=skiroc.TOARise(ichan);
      m_toa_fall[ichan]=skiroc.TOAFall(ichan);
    }
    m_tree->Fill();
  }
}

void TreeProducer::beginJob()
{
}

void TreeProducer::endJob()
{
}

void TreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TreeProducer);
