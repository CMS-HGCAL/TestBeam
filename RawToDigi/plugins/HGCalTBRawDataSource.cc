#include "HGCal/RawToDigi/plugins/HGCalTBRawDataSource.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

HGCalTBRawDataSource::HGCalTBRawDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :
  edm::ProducerSourceFromFiles(pset, desc, true),
  m_electronicMap(pset.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt")),
  m_outputCollectionName(pset.getUntrackedParameter<std::string>("OutputCollectionName","SKIROC2CMSDATA")),
  m_nOrmBoards(pset.getUntrackedParameter<unsigned int>("NOrmBoards", 1)),
  m_nHexaboards(pset.getUntrackedParameter<unsigned int>("NHexaBoards", 1)),
  m_nWords(pset.getUntrackedParameter<unsigned int> ("NumberOf32BitsWordsPerReadOut",30788))
{
  produces<HGCalTBSkiroc2CMSCollection>(m_outputCollectionName);
  
  m_event = 0;
  m_fileId=0;
 
  m_skiMask = 0xffffffff;
  switch( m_nHexaboards ){
  case 1 : m_skiMask = 0x0f000000;break;
  case 2 : m_skiMask = 0xff000000;break;
  case 3 : m_skiMask = 0xff0f0000;break;
  case 4 : m_skiMask = 0xffff0000;break;
  case 5 : m_skiMask = 0xffff0f00;break;
  case 6 : m_skiMask = 0xffffff00;break;
  case 7 : m_skiMask = 0xffffff0f;break;
  case 8 : m_skiMask = 0xffffffff;break;
  default : { std::cout << "wrong skiroc number: " << m_nHexaboards << std::endl; exit(1); }
  }

  m_buffer=new char[m_nWords*4];
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }

  std::cout << pset << std::endl;

}

bool HGCalTBRawDataSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{	
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
  }
  m_input.seekg( m_event*m_nWords*4, std::ios::beg );
  m_input.read ( m_buffer, m_nWords*4 );
  if( !m_input.good() ){ m_input.close(); m_fileId++; m_event=0; return setRunAndEventInfo(id, time, evType);}
  std::vector<uint32_t> rawData;
  for( size_t i=0; i<m_nWords; i++ ){
    uint32_t aint;
    char buf[] = {m_buffer[i*4+3],m_buffer[i*4+2],m_buffer[i*4+1],m_buffer[i*4]};
    memcpy(&aint, &buf, sizeof(aint));
    rawData.push_back(aint);
  }

  m_decodedData=decode_raw_32bit(rawData);
  return true;
}

std::vector< std::array<uint16_t,1924> > HGCalTBRawDataSource::decode_raw_32bit(std::vector<uint32_t> &raw){
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
      int skiId;
      if( iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA%2==0 )//not flipped
	skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*(iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iski)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
      else
	skiId = HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*(iski/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-iski%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA);
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
    skirocs->push_back(skiroc);
  }
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
  desc.addUntracked<unsigned int>("NOrmBoards");
  desc.addUntracked<unsigned int>("NHexaBoards");
  desc.addUntracked<unsigned int>("NumberOf32BitsWordsPerReadOut");
  descriptions.add("source", desc);
}

void HGCalTBRawDataSource::endJob()
{
  //delete m_buffer;
}

#include "FWCore/Framework/interface/InputSourceMacros.h"
DEFINE_FWK_INPUT_SOURCE(HGCalTBRawDataSource);
