#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
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
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include <iomanip>
#include <set>

//#define DEBUG

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
  bool m_writeNoisyChannelFile;
  std::string m_noisyChannels_filename;
  bool m_writeTreeOutput;
  std::string m_electronicMap;
  int m_NTSForPedestalComputation;

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
  m_writeNoisyChannelFile(iConfig.getUntrackedParameter<bool>("WriteNoisyChannelsFile",false)),
  m_noisyChannels_filename( iConfig.getUntrackedParameter<std::string>("NoisyChannelsFileName",std::string("noisyChannels.txt")) ),
  m_writeTreeOutput(iConfig.getUntrackedParameter<bool>("WriteTreeOutput",false)),
  m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt")),
  m_NTSForPedestalComputation(iConfig.getUntrackedParameter<int>("NTSForPedestalComputation",1))    
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
    //std::cout << skiroc << std::endl;
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
	if( rollpositions[it]<m_NTSForPedestalComputation ){ //consider only a certain number of time samples for pedestal subtraction
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
    h->SetOption("colztext");
    InitTH2Poly(*h, ib, 0, 0);
    chanMap.insert( std::pair<int,TH2Poly*>(ib,h) );
    TFileDirectory hgpdir = dir.mkdir( "HighGainPedestal" );
    TFileDirectory lgpdir = dir.mkdir( "LowGainPedestal" );
    TFileDirectory hgndir = dir.mkdir( "HighGainNoise" );
    TFileDirectory lgndir = dir.mkdir( "LowGainNoise" );
    for( size_t it=0; it<NUMBER_OF_SCA; it++ ){

      h=hgpdir.make<TH2Poly>();
      os.str("");
      os<<"SCA_"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, ib, 0, 0);
      hgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgpdir.make<TH2Poly>();
      os.str("");
	
      os<<"SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, ib, 0, 0);
      lgMeanMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=hgndir.make<TH2Poly>();
      os.str("");
      os<<"SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, ib, 0, 0);
      hgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );

      h=lgndir.make<TH2Poly>();
      os.str("");
      os<<"SCA"<<it;
      h->SetName(os.str().c_str());
      h->SetTitle(os.str().c_str());
      h->SetOption("colztext");
      InitTH2Poly(*h, ib, 0, 0);
      lgRMSMap.insert( std::pair<int,TH2Poly*>(100*ib+it,h) );
    }
  }
  
  for( std::map<int,hgcal_channel>::iterator it=m_channelMap.begin(); it!=m_channelMap.end(); ++it ){
    std::sort( it->second.highGain.begin(),it->second.highGain.end() );
    std::sort( it->second.lowGain.begin(),it->second.lowGain.end() );
    unsigned int size = it->second.highGain.size();
    int medianIndex = int( 0.5*(size-1) );
    it->second.medianHG = it->second.highGain.at(medianIndex) ;
    it->second.medianLG = it->second.lowGain.at(medianIndex) ;
    int sigma_1_Index = int( 0.16*(size-1) );
    int sigma_3_Index = int( 0.84*(size-1) );
    it->second.rmsHG=0.5*(it->second.highGain.at(sigma_3_Index)-it->second.highGain.at(sigma_1_Index));
    it->second.rmsLG=0.5*(it->second.lowGain.at(sigma_3_Index)-it->second.lowGain.at(sigma_1_Index));    
  }

  
  std::fstream pedestalHG;
  std::fstream pedestalLG;
  if( m_writePedestalFile ){
    pedestalHG.open(m_pedestalHigh_filename,std::ios::out);
    pedestalLG.open(m_pedestalLow_filename,std::ios::out);
  }


  TTree* pedestalTree = NULL;
  int iboard, iski, ichan, sca;
  double median_high, RMS_high, median_low, RMS_low;
  if( m_writeTreeOutput ){
    pedestalTree=fs->make<TTree>("pedestalTree", "pedestalTree");
    pedestalTree->Branch("board", &iboard);
    pedestalTree->Branch("skiroc", &iski);
    pedestalTree->Branch("channel", &ichan);
    pedestalTree->Branch("sca", &sca);
    pedestalTree->Branch("median_high", &median_high);
    pedestalTree->Branch("RMS_high", &RMS_high);
    pedestalTree->Branch("median_low", &median_low);
    pedestalTree->Branch("RMS_low", &RMS_low);
  }

  std::fstream noisyChannels;
  if( m_writeNoisyChannelFile )
    noisyChannels.open(m_noisyChannels_filename,std::ios::out);
  std::map<int,double> meanNoise;
  std::map<int,double> rmsNoise;
  std::map<int,int> countNoise;

  for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
    iboard=(*it).first/1000;
    iski=((*it).first%1000)/100;
    ichan=(*it).first%100;
    if( m_writePedestalFile ){
      pedestalHG << iboard << " " << iski << " " << ichan ;
      pedestalLG << iboard << " " << iski << " " << ichan ;
    }
    HGCalTBDetId detid=(*it).second;
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( detid.layer(), 0, 0, detid.iu(), detid.iv(), m_sensorsize );


    double iux = CellCentreXY.first;
    double iuy = CellCentreXY.second;

    
    for( size_t it=0; it<NUMBER_OF_SCA; it++ ){

      int key=iboard*100000+iski*10000+ichan*100+it;
      std::map<int,hgcal_channel>::iterator iter=m_channelMap.find(key);
      float hgMean=iter->second.medianHG;
      float lgMean=iter->second.medianLG;
      float hgRMS=iter->second.rmsHG;
      float lgRMS=iter->second.rmsLG;

      hgMeanMap[ 100*iboard+it ]->Fill(iux , iuy, hgMean );
      lgMeanMap[ 100*iboard+it ]->Fill(iux , iuy, lgMean );
      hgRMSMap[ 100*iboard+it ]->Fill(iux , iuy, hgRMS );
      lgRMSMap[ 100*iboard+it ]->Fill(iux , iuy, lgRMS );
      
      if( m_writeTreeOutput ){
	sca = it;
	median_high = hgMean;
	RMS_high = hgRMS;
	median_low = lgMean;
	RMS_low = lgRMS;
	pedestalTree->Fill();
      }

      if( m_writePedestalFile ){
      	pedestalHG << " " << hgMean << " " << hgRMS;
	pedestalLG << " " << lgMean << " " << lgRMS;
      }
    }
    chanMap[ iboard ]->Fill(iux , iuy, iski*1000+ichan );
    if( m_writePedestalFile ){
      pedestalHG << std::endl;
      pedestalLG << std::endl;
    }
    if( m_writeNoisyChannelFile ){
      if( detid.cellType()!=0 && detid.cellType()!=5  && detid.cellType()!=4 )continue;
      int key=iboard*100000+iski*10000+ichan*100;
      std::map<int,hgcal_channel>::iterator iter=m_channelMap.find(key);
      float hgRMS=iter->second.rmsHG;
      if( meanNoise.find(iboard)==meanNoise.end() ){
	meanNoise[iboard]=hgRMS;
	rmsNoise[iboard]=hgRMS*hgRMS;
	countNoise[iboard]=1;
      }
      else{
	meanNoise[iboard]+=hgRMS;
	rmsNoise[iboard]+=hgRMS*hgRMS;
	countNoise[iboard]+=1;
      }
    }
  }
  if( m_writePedestalFile ){
    pedestalHG.close();
    pedestalLG.close();
  }
  if( m_writeNoisyChannelFile ){
    for( std::map<int,double>::iterator it=meanNoise.begin(); it!=meanNoise.end(); ++it ){
      meanNoise[ it->first ] = meanNoise[ it->first ]/countNoise[ it->first ];
      rmsNoise[ it->first ] = std::sqrt( rmsNoise[ it->first ]/countNoise[ it->first ] - meanNoise[ it->first ]*meanNoise[ it->first ] );
      std::cout << it->first << " " << meanNoise[ it->first ] << " " << rmsNoise[ it->first ] << std::endl;
    }
    for( std::set< std::pair<int,HGCalTBDetId> >::iterator it=setOfConnectedDetId.begin(); it!=setOfConnectedDetId.end(); ++it ){
      HGCalTBDetId detid=(*it).second;
      if( detid.cellType()!=0 && detid.cellType()!=4 )continue;
      int iboard=(*it).first/1000;
      int iski=((*it).first%1000)/100;
      int ichan=(*it).first%100;
      int key=iboard*100000+iski*10000+ichan*100;//we use SCA 0
      std::map<int,hgcal_channel>::iterator iter=m_channelMap.find(key);
      if( iter->second.rmsHG-meanNoise[iboard]>3*rmsNoise[iboard] )
	noisyChannels << iboard << " " << iski << " " << ichan << std::endl;
    }
    noisyChannels.close();
  }
}

void PedestalPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PedestalPlotter);
