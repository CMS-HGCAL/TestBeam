#include "HGCal/RawToDigi/plugins/HGCalTBRawDataSource.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

HGCalTBRawDataSource::HGCalTBRawDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :
  edm::ProducerSourceFromFiles(pset, desc, true),
  m_filePath(pset.getUntrackedParameter<std::string>("FilePath", "./data/")),
  m_fileName(pset.getUntrackedParameter<std::string>("FileName", "HexaData_Run0000.raw")),
  m_electronicMap(pset.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt")),
  m_outputCollectionName(pset.getUntrackedParameter<std::string>("OutputCollectionName","SKIROC2CMSDATA")),
  m_nOrmBoards(pset.getUntrackedParameter<unsigned int>("NOrmBoards", 1)),
  m_nHexaboards(pset.getUntrackedParameter<unsigned int>("NHexaBoards", 1)),
  m_nSkirocsPerHexa(pset.getUntrackedParameter<unsigned int>("NSkirocsPerHexa", 4)),
  m_nChannelsPerSkiroc(pset.getUntrackedParameter<unsigned int>("NChannelsPerSkiroc", 64)),
  m_nWords(pset.getUntrackedParameter<unsigned int> ("NumberOf32BitsWordsPerReadOut",30788))
{
  std::cout << "Hey je commence" << std::endl;
  sleep(10);

  produces<HGCalTBSkiroc2CMSCollection>();
  //produces<RunData>("RunData");
  
  m_event = 0;
  m_fileId=0;
 
  m_skiMask = 0xffffffff;
  switch( m_nOrmBoards ){
  case 1 : m_skiMask = 0x0f000000;break;
  case 2 : m_skiMask = 0xff000000;break;
  case 3 : m_skiMask = 0xff0f0000;break;
  case 4 : m_skiMask = 0xffff0000;break;
  case 5 : m_skiMask = 0xffff0f00;break;
  case 6 : m_skiMask = 0xffffff00;break;
  case 7 : m_skiMask = 0xffffff0f;break;
  case 8 : m_skiMask = 0xffffffff;break;
  default : { std::cout << "wrong skiroc number: " << m_nOrmBoards << std::endl; exit(1); }
  }

  m_buffer=new char[m_nWords*4];
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }

}

bool HGCalTBRawDataSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{	

  if( m_event==0 ){
    m_input.open(m_fileName.c_str(), std::ios::in|std::ios::binary);
    if( !m_input.is_open() ){
      std::cout << "PROBLEM : Not able to open " << m_fileName << "\t -> return false (so end of process I guess)" << std::endl;
      return false;
    }
  }
  if( !m_input.good() ) return false;
  std::cout << "start reading event " << m_event << std::endl;
  m_input.read ( m_buffer, m_nWords*4 );
  m_input.seekg( (m_event+1)*m_nWords*4, std::ios::beg );
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
  std::cout<<"In decoder"<<std::endl;
  std::cout<<"SkiMask = " << std::setw(8) << std::setfill('0') << std::hex << m_skiMask << "\t Length of Raw " << std::dec << raw.size() << std::endl;;
    // Check that an external mask agrees with first 32-bit word in data
  if (m_skiMask!=raw[0])
    std::cout << "You extarnal mask (" << std::hex << raw[0] << ") does not agree with the one found in data (" << std::hex << m_skiMask << ")" << std::endl;
  
  for (int b=0; b<20; b++){
    std::cout << std::dec << b << "\t";
    std::cout << std::setw(8) << std::setfill('0') << std::hex << raw[b] <<std::endl;
  }
  std::cout<< std::dec << 30786 << "\t" << std::setw(8) << std::setfill('0') << std::hex << raw[30786] <<std::endl;
  std::cout<< std::dec << 30787 << "\t" << std::setw(8) << std::setfill('0') << std::hex << raw[30787] <<std::endl;

  // First, we need to determine how many skiRoc data is presnt
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
    std::array<HGCalTBDetId,64> detids;
    for (size_t ichan = 0; ichan < m_nChannelsPerSkiroc; ichan++) {
      HGCalTBElectronicsId eid(iski, ichan);
      if (!essource_.emap_.existsEId(eid.rawId())) {
	std::cout << eid.rawId() << " is not a correct id, or can not be found in electronics map" << std::endl;
	//exit(1);
      }
      else{ 
	HGCalTBDetId did = essource_.emap_.eid2detId(eid);
	detids[ichan]=did;
      }
    }
    HGCalTBSkiroc2CMS skiroc( m_decodedData.at(iski),detids);
    if(!skiroc.check())
      exit(1);
    std::cout << skiroc << std::endl;
    getchar();
    skirocs->push_back(skiroc);
  }
  e.put(skirocs, m_outputCollectionName);
  
  m_event++;
  
}

void HGCalTBRawDataSource::endJob()
{
  delete m_buffer;
}
