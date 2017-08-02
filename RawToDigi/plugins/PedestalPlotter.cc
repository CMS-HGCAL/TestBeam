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

struct hgcal_channel{
  hgcal_channel() : key(0),
		    medianHG(0.),
		    medianLG(0.),
		    rmsHG(0.),
		    rmsLG(0.){;}
  int key;
  float medianHG;
  float medianLG;
  float rmsHG;
  float rmsLG;
  std::vector<float> highGain;
  std::vector<float> lowGain;
};

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
  std::map<int,hgcal_channel> m_channelMap;

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
  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);

  m_numberOfBoards = skirocs->size()/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA;

  m_evtID++;

  if( !skirocs->size() ) return;
  
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
	if( rollpositions[it]<2 ){ //keep only 2 first time sample for pedestal evaluation
	  uint32_t key=iboard*100000+(iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)*10000+ichan*100+it;
	  std::map<int,hgcal_channel>::iterator iter=m_channelMap.find(key);
	  if( iter==m_channelMap.end() ){
	    hgcal_channel tmp;
	    tmp.key=key;
	    std::vector<float> vecH,vecL;
	    vecH.push_back(skiroc.ADCHigh(ichan,it));
	    vecL.push_back(skiroc.ADCLow(ichan,it));
	    tmp.highGain=vecH;
	    tmp.lowGain=vecL;
	    std::pair<int,hgcal_channel> p(key,tmp);
	    m_channelMap.insert( p );
	  }
	  else{
	    iter->second.highGain.push_back(skiroc.ADCHigh(ichan,it));
	    iter->second.lowGain.push_back(skiroc.ADCLow(ichan,it));
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
  std::map<int,TH2Poly*>  chanMap;
  std::ostringstream os( std::ostringstream::ate );
  TH2Poly *h;
  for(size_t ib = 0; ib<m_numberOfBoards; ib++) {
    os.str("");
    os << "HexaBoard" << ib ;
    TFileDirectory dir = fs->mkdir( os.str().c_str() );
    h=dir.make<TH2Poly>();
    os.str("");
    os<<"ChannelMapping";
    h->SetName(os.str().c_str());
    h->SetTitle(os.str().c_str());
    InitTH2Poly(*h, (int)ib, 0, 0);
    chanMap.insert( std::pair<int,TH2Poly*>(ib,h) );
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
  
  for( std::map<int,hgcal_channel>::iterator it=m_channelMap.begin(); it!=m_channelMap.end(); ++it ){
    std::sort( it->second.highGain.begin(),it->second.highGain.end() );
    std::sort( it->second.lowGain.begin(),it->second.lowGain.end() );
    unsigned int size = it->second.highGain.size();
    int medianIndex = size%2==0 ? size/2-1 : size/2 ;
    it->second.medianHG = it->second.highGain.at(medianIndex) ;
    it->second.medianLG = it->second.lowGain.at(medianIndex) ;
    for( unsigned int ii=medianIndex; ii<size; ii++ ){
      it->second.rmsHG+=(it->second.medianHG-it->second.highGain.at(ii))*(it->second.medianHG-it->second.highGain.at(ii));
      it->second.rmsLG+=(it->second.medianLG-it->second.lowGain.at(ii))*(it->second.medianLG-it->second.lowGain.at(ii));
    }
    it->second.rmsHG=std::sqrt(2*it->second.rmsHG/size);
    it->second.rmsLG=std::sqrt(2*it->second.rmsLG/size);
  }

  
  std::fstream pedestalHG;
  std::fstream pedestalLG;
  if( m_writePedestalFile ){
    pedestalHG.open(m_pedestalHigh_filename,std::ios::out);
    pedestalLG.open(m_pedestalLow_filename,std::ios::out);
  }
  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    int iboard=(*it).first/1000;
    int iski=((*it).first%1000)/100;
    int ichan=(*it).first%100;
    if( m_writePedestalFile ){
      pedestalHG << iboard << " " << iski << " " << ichan ;
      pedestalLG << iboard << " " << iski << " " << ichan ;
    }
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), detid.sensorIU(), detid.sensorIV(), detid.iu(), detid.iv(), m_sensorsize );
    double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.first - HGCAL_TB_GEOMETRY::DELTA) ;
    double iuy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + HGCAL_TB_GEOMETRY::DELTA) : (CellCentreXY.second - HGCAL_TB_GEOMETRY::DELTA);
    for( size_t it=0; it<NUMBER_OF_SCA; it++ ){
      int key=iboard*100000+iski*10000+ichan*100+it;
      std::map<int,hgcal_channel>::iterator iter=m_channelMap.find(key);
      float hgMean=iter->second.medianHG;
      float lgMean=iter->second.medianLG;
      float hgRMS=iter->second.rmsHG;
      float lgRMS=iter->second.rmsLG;
      hgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgMean );
      lgMeanMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgMean );
      hgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, hgRMS );
      lgRMSMap[ 100*iboard+it ]->Fill(iux/2 , iuy, lgRMS );
      if( m_writePedestalFile ){
	pedestalHG << " " << hgMean << " " << hgRMS;
	pedestalLG << " " << lgMean << " " << lgRMS;;
      }
    }
    chanMap[ iboard ]->Fill(iux/2 , iuy, iski*1000+ichan );
    if( m_writePedestalFile ){
      pedestalHG << std::endl;
      pedestalLG << std::endl;
    }
  }
  if( m_writePedestalFile ){
    pedestalHG.close();
    pedestalLG.close();
  }
}

void PedestalPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PedestalPlotter);
