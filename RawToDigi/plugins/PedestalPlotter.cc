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
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include <iomanip>
#include <set>

class PedestalPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit PedestalPlotter(const edm::ParameterSet&);
  ~PedestalPlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  int m_sensorsize;

  bool m_writePedestalFile;
  std::string m_pedestalHigh_filename;
  std::string m_pedestalLow_filename;
  std::string m_electronicMap;

  int m_evtID;
  uint16_t m_numberOfBoards;
  std::map<int,float> m_meanHGMap;
  std::map<int,float> m_meanLGMap;
  std::map<int,float> m_rmsHGMap;
  std::map<int,float> m_rmsLGMap;
  std::map<int,int> m_counterMap;
  // std::map<int,TH1F*> m_h_adcHigh;
  // std::map<int,TH1F*> m_h_adcLow;

  edm::EDGetTokenT<HGCalTBSkiroc2CMSCollection> m_HGCalTBSkiroc2CMSCollection;

  HGCalTBTopology IsCellValid;
  HGCalTBCellVertices TheCell;
  std::vector<std::pair<double, double>> CellXY;
  std::pair<double, double> CellCentreXY;
  std::set< std::pair<int,HGCalTBDetId> > setOfConnectedDetId;
};

PedestalPlotter::PedestalPlotter(const edm::ParameterSet& iConfig) :
  m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
  m_writePedestalFile(iConfig.getUntrackedParameter<bool>("WritePedestalFile",false)),
  m_pedestalHigh_filename( iConfig.getUntrackedParameter<std::string>("HighGainPedestalFileName",std::string("pedestalHG.txt")) ),
  m_pedestalLow_filename( iConfig.getUntrackedParameter<std::string>("LowGainPedestalFileName",std::string("pedestalLG.txt")) ),
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt"))
{
  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

  m_evtID=0;

  // usesResource("TFileService");
  // edm::Service<TFileService> fs;
  // std::ostringstream os( std::ostringstream::ate );
  // TH1F* htmp1;
  // for(size_t ib = 0; ib<HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
  //   std::cout << "Hexaboard " << ib << std::endl;
  //   for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
  //     os.str("");os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
  //     TFileDirectory dir = fs->mkdir( os.str().c_str() );
  //     for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
  // 	for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
  // 	  os.str("");
  // 	  os << "HighGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_SCA" << it ;
  // 	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
  // 	  m_h_adcHigh.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
  // 	  os.str("");
  // 	  os << "LowGain_HexaBoard" << ib << "_Chip" << iski << "_Channel" << ichan << "_SCA" << it ;
  // 	  htmp1=dir.make<TH1F>(os.str().c_str(),os.str().c_str(),1000,0,4096);
  // 	  m_h_adcLow.insert( std::pair<int,TH1F*>(ib*100000+iski*10000+ichan*100+it, htmp1) );
  // 	}
  //     }
  //   }
  // }
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  };
  std::cout << iConfig.dump() << std::endl;
}


PedestalPlotter::~PedestalPlotter()
{

}

void PedestalPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  //  std::cout << "Here I am" << std::endl;
  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);

  //std::cout << skirocs->size() << std::endl;

  m_numberOfBoards = skirocs->size()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;

  // std::cout << "m_numberOfBoards = " << m_numberOfBoards << std::endl;

  m_evtID++;

  for( size_t iski=0;iski<skirocs->size(); iski++ ){
    HGCalTBSkiroc2CMS skiroc=skirocs->at(iski);
    std::vector<int> rollpositions=skiroc.rollPositions();
    int iboard=iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;
    for( size_t ichan=0; ichan<HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
      HGCalTBDetId detid=skiroc.detid( ichan );
      HGCalTBElectronicsId eid( essource_.emap_.detId2eid(detid.rawId()) );
      if( essource_.emap_.existsEId(eid) ){
	std::pair<int,HGCalTBDetId> p( iboard*1000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*100+ichan,skiroc.detid(ichan) );
	setOfConnectedDetId.insert(p);
      }
      else continue;
      for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
	if( rollpositions[it]<9 ){ //rm on track samples+2 last time sample which show weird behaviour
	  uint32_t key=iboard*100000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*10000+ichan*100+it;
	  // m_h_adcHigh[key]->Fill(skiroc.ADCHigh(ichan,it));
	  // m_h_adcLow[key]->Fill(skiroc.ADCLow(ichan,it));
	  if( m_meanHGMap.find(key)==m_meanHGMap.end() ){
	    m_meanHGMap[key]=skiroc.ADCHigh(ichan,it);
	    m_meanLGMap[key]=skiroc.ADCLow(ichan,it);
	    m_rmsHGMap[key]=skiroc.ADCHigh(ichan,it)*skiroc.ADCHigh(ichan,it);
	    m_rmsLGMap[key]=skiroc.ADCLow(ichan,it)*skiroc.ADCLow(ichan,it);
	    m_counterMap[key]=1;
	  }
	  else{
	    m_meanHGMap[key]+=skiroc.ADCHigh(ichan,it);
	    m_meanLGMap[key]+=skiroc.ADCLow(ichan,it);
	    m_rmsHGMap[key]+=skiroc.ADCHigh(ichan,it)*skiroc.ADCHigh(ichan,it);
	    m_rmsLGMap[key]+=skiroc.ADCLow(ichan,it)*skiroc.ADCLow(ichan,it);
	    m_counterMap[key]+=1;
	  }
	}
      }

    }
  }
}

void PedestalPlotter::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
  double HexY[HGCAL_TB_GEOMETRY::MAXVERTICES] = {0.};
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

void PedestalPlotter::beginJob()
{
}

void PedestalPlotter::endJob()
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  std::map<int,TH2Poly*>  hgMeanMap;
  std::map<int,TH2Poly*>  lgMeanMap;
  std::map<int,TH2Poly*>  hgRMSMap;
  std::map<int,TH2Poly*>  lgRMSMap;
  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h;
  for(size_t ib = 0; ib<m_numberOfBoards; ib++) {
    os.str("");
    os << "HexaBoard" << ib ;
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    TFileDirectory hgpdir = dir.mkdir( "HighGainPedestal" );
    TFileDirectory lgpdir = dir.mkdir( "LowGainPedestal" );
    TFileDirectory hgndir = dir.mkdir( "HighGainNoise" );
    TFileDirectory lgndir = dir.mkdir( "LowGainNoise" );
    for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
      h=hgpdir.make<TH2Poly>();
      os.str("");
      os<<"SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      hgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgpdir.make<TH2Poly>();
      os.str("");
      os<<"SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      lgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=hgndir.make<TH2Poly>();
      os.str("");
      os<<"SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      hgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgndir.make<TH2Poly>();
      os.str("");
      os<<"SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      InitTH2Poly(*h, (int)ib, 0, 0);
      lgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
    }
  }

  std::fstream pedestalHG;pedestalHG.open(m_pedestalHigh_filename,std::ios::out);
  std::fstream pedestalLG;pedestalLG.open(m_pedestalLow_filename,std::ios::out);
  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    pedestalHG << iboard << " " << iski << " " << ichan ;
    pedestalLG << iboard << " " << iski << " " << ichan ;
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA) ;
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
    for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
      int key=iboard*100000+iski*10000+ichan*100+it;
      float hgMean=m_meanHGMap[key]/m_counterMap[key];
      float lgMean=m_meanLGMap[key]/m_counterMap[key];
      float hgRMS=std::sqrt(m_rmsHGMap[key]/m_counterMap[key]-m_meanHGMap[key]/m_counterMap[key]*m_meanHGMap[key]/m_counterMap[key]);
      float lgRMS=std::sqrt(m_rmsLGMap[key]/m_counterMap[key]-m_meanLGMap[key]/m_counterMap[key]*m_meanLGMap[key]/m_counterMap[key]);
      hgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgMean );
      lgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgMean );
      hgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgRMS );
      lgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgRMS );
      pedestalHG << " " << hgMean << " " << hgRMS;
      pedestalLG << " " << lgMean << " " << lgRMS;;
    }
  }
  pedestalHG.close();
  pedestalLG.close();
}

void PedestalPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PedestalPlotter);
