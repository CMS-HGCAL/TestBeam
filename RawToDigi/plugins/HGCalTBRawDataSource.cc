#include "HGCal/RawToDigi/plugins/HGCalTBRawDataSource.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include <iomanip>
#include <ctime>
#include <cmath>
HGCalTBRawDataSource::HGCalTBRawDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :
  edm::ProducerSourceFromFiles(pset, desc, true),
  m_electronicMap(pset.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt")),
  m_outputCollectionName(pset.getUntrackedParameter<std::string>("OutputCollectionName","SKIROC2CMSDATA")),
  m_nWords(pset.getUntrackedParameter<unsigned int> ("NumberOf32BitsWordsPerReadOut",30788)),
  m_headerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheHeader",8)),
  m_trailerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheTrailer",4)),
  m_eventTrailerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheEventTrailers",4))
{
  produces<HGCalTBSkiroc2CMSCollection>(m_outputCollectionName);
  
  m_event = 0;
  m_fileId=0;
  m_meanReadingTime=0;
  m_rmsReadingTime=0;

  m_buffer=new char[m_nWords*4+m_eventTrailerSize];//4 bytes per 32 bits
  m_header=new char[m_headerSize];
  m_trailer=new char[m_trailerSize];
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }
  
  std::cout << pset << std::endl;

}

bool HGCalTBRawDataSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{
  clock_t start;
  start=clock();
  if( m_fileId == fileNames().size() ) return false;
  if (fileNames()[m_fileId] != "file:DUMMY")
    m_fileName = fileNames()[m_fileId];
  if (m_fileName.find("file:") == 0) m_fileName = m_fileName.substr(5);
  if( m_event==0 ){
    m_input.open(m_fileName.c_str(), std::ios::in|std::ios::binary);
    if( !m_input.is_open() ){
      std::cout << "PROBLEM : Not able to open " << m_fileName << "\t -> return false (so end of process I guess)" << std::endl;
      return false;
    }
    m_input.seekg( 0, std::ios::beg );
    m_input.read ( m_header, m_headerSize );
    char buf0[] = {m_header[0],m_header[1],m_header[2],m_header[3]};
    memcpy(&timeStartRun, &buf0, sizeof(timeStartRun));
    char buf1[] = {m_header[4],m_header[5],m_header[6],m_header[7]};
    uint32_t aint;
    memcpy(&aint, &buf1, sizeof(aint));
    m_nOrmBoards=aint & 0xff;
    m_run=aint >> 8 ;
  }

  m_decodedData.clear();
  std::vector<uint32_t> rawData;
  for( uint16_t iorm=0; iorm<m_nOrmBoards; iorm++ ){
    m_input.seekg( m_headerSize+m_nOrmBoards*m_event*(m_nWords*4+m_eventTrailerSize)+iorm*(m_nWords*4+m_eventTrailerSize), std::ios::beg );
    m_input.read ( m_buffer, m_nWords*4+m_eventTrailerSize );
    if( !m_input.good() ){
      m_input.close();
      m_fileId++;
      if( (uint32_t)(m_fileId+1)<fileNames().size() )
	m_event=0;
      return setRunAndEventInfo(id, time, evType);
    }
    for( size_t i=0; i<m_nWords; i++ ){
      uint32_t aint;
      char buf[] = {m_buffer[i*4+3],m_buffer[i*4+2],m_buffer[i*4+1],m_buffer[i*4]};
      memcpy(&aint, &buf, sizeof(aint));
      rawData.push_back(aint);
    }
    uint32_t evtTrailer;
    char buf[] = {m_buffer[m_nWords*4],m_buffer[m_nWords*4+1],m_buffer[m_nWords*4+2],m_buffer[m_nWords*4+3]};
    memcpy(&evtTrailer, &buf, sizeof(evtTrailer));
    uint32_t ormId=evtTrailer&0xff;
    uint32_t evtNumber=evtTrailer>>0x8;
    if( ormId != iorm )
      std::cout << "Problem in event trailer : wrong ORM id -> evtTrailer&0xff = " << std::dec << ormId << "\t iorm = " << iorm << std::endl;
    if( evtNumber != m_event+1 )
      std::cout << "Problem in event trailer : wrong event number -> evtTrailer>>8 = " << std::dec << evtNumber << "\t m_event+1 = " << m_event+1 << std::endl;
    std::vector< std::array<uint16_t,1924> > decodedData=decode_raw_32bit(rawData);
    m_decodedData.insert(m_decodedData.end(),decodedData.begin(),decodedData.end());
  }
  clock_t end;
  end=clock();
  m_meanReadingTime += (float)(end-start)/CLOCKS_PER_SEC;
  m_rmsReadingTime += (float)(end-start)/CLOCKS_PER_SEC*(end-start)/CLOCKS_PER_SEC;
  return true;
}

std::vector< std::array<uint16_t,1924> > HGCalTBRawDataSource::decode_raw_32bit(std::vector<uint32_t> &raw){
  m_skiMask=raw[0];
  const std::bitset<32> ski_mask(m_skiMask);
  const int mask_count = ski_mask.count();
  std::vector< std::array<uint16_t, 1924> > skirocs(mask_count, std::array<uint16_t, 1924>());
  uint32_t x;
  const int offset = 1; // Due to FF or other things in data head
  for(int  i = 0; i < 1924; i++){
    for (int j = 0; j < 16; j++){
      x = raw[offset + i*16 + j];
      int k = 0;
      for (int fifo = 0; fifo < 32; fifo++){
	if (m_skiMask & (1<<fifo)){
	  skirocs[k][i] = skirocs[k][i] | (uint16_t) (((x >> fifo ) & 1) << (15 - j));
	  k++;
	}
      }
    }
  }
  return skirocs;
}

void HGCalTBRawDataSource::produce(edm::Event & e)
{
  std::auto_ptr<HGCalTBSkiroc2CMSCollection> skirocs(new HGCalTBSkiroc2CMSCollection);

  for( size_t iski=0; iski<m_decodedData.size(); iski++){
    std::vector<HGCalTBDetId> detids;
    for (size_t ichan = 0; ichan < HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++) {
      int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*(iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iski)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
      HGCalTBElectronicsId eid(skiId,ichan);      
      if (!essource_.emap_.existsEId(eid.rawId())) {
	HGCalTBDetId did(-1);
	detids.push_back(did);
      }
      else{ 
	HGCalTBDetId did = essource_.emap_.eid2detId(eid);
	detids.push_back(did);
      }
    }
    std::vector<uint16_t> vdata;vdata.insert(vdata.end(),m_decodedData.at(iski).begin(),m_decodedData.at(iski).end());
    HGCalTBSkiroc2CMS skiroc( vdata,detids );
    if(!skiroc.check())
      exit(1);
    //std::cout << skiroc << std::endl;
    skirocs->push_back(skiroc);
  }
  //std::cout << "skirocs->size() = " << skirocs->size() << std::endl;
  e.put(skirocs, m_outputCollectionName);
  m_event++;
}

void HGCalTBRawDataSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setComment("TEST");
  desc.addUntracked<int>("run", 101);
  desc.addUntracked<std::vector<std::string> >("fileNames");
  desc.addUntracked<std::string>("ElectronicMap");
  desc.addUntracked<std::string>("OutputCollectionName");
  desc.addUntracked<unsigned int>("NumberOf32BitsWordsPerReadOut");
  desc.addUntracked<unsigned int> ("NumberOfBytesForTheHeader");
  desc.addUntracked<unsigned int> ("NumberOfBytesForTheTrailer");
  desc.addUntracked<unsigned int> ("NumberOfBytesForTheEventTrailers");
  descriptions.add("source", desc);
}

void HGCalTBRawDataSource::endJob()
{
  std::cout << "mean reading = " << m_meanReadingTime/m_event
	    << "  +-  " << std::sqrt(m_rmsReadingTime/m_event-m_meanReadingTime/m_event*m_meanReadingTime/m_event)
	    << std::endl;

  //delete m_buffer;
}

#include "FWCore/Framework/interface/InputSourceMacros.h"
DEFINE_FWK_INPUT_SOURCE(HGCalTBRawDataSource);
