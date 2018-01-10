#include "HGCal/RawToDigi/plugins/HGCalTBRawDataSource.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include <stdlib.h>
#include <iomanip>
#include <ctime>
#include <cmath>

HGCalTBRawDataSource::HGCalTBRawDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :
  edm::ProducerSourceFromFiles(pset, desc, true),
  m_electronicMap(pset.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt")),
  m_outputCollectionName(pset.getUntrackedParameter<std::string>("OutputCollectionName","SKIROC2CMSDATA")),
  m_nWords(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesPerReadOut",30784)),
  m_headerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheHeader",48)),
  m_eventTrailerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheEventTrailers",2)),
  m_nSkipEvents(pset.getUntrackedParameter<unsigned int> ("NSkipEvents",1)),
  m_compressedData(pset.getUntrackedParameter<bool> ("CompressedData",false))
{
  produces<HGCalTBSkiroc2CMSCollection>(m_outputCollectionName);
  
  m_event = 0;
  m_fileId=0;

  m_buffer=new char[m_nWords+m_eventTrailerSize];
  m_header=new char[m_headerSize];
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
    m_input.seekg( 0, std::ios::beg );
    m_input.read ( m_header, m_headerSize );
    std::vector<uint16_t> bitstring;
    for( int i=0; i<int(m_headerSize/2); i++ ){
      uint16_t aint=0;
      char buf[] = {m_header[i*2],m_header[i*2+1]};
      memcpy(&aint, &buf, sizeof(aint));
      bitstring.push_back(aint);
    }
    std::cout << "bit string = \t";
    for( std::vector<uint16_t>::iterator it=bitstring.begin(); it!= bitstring.end(); ++it )
      std::cout << " " << std::setfill('0') << std::setw(4) << std::hex << (*it) ;
    std::cout << std::endl;
  }

  uint64_t nBytesToSkip=m_headerSize+m_nSkipEvents*(m_nWords+m_eventTrailerSize);
  std::vector<uint16_t> rawData;
  m_input.seekg( (std::streamoff)nBytesToSkip+(std::streamoff)m_event*(m_nWords+m_eventTrailerSize) );
  m_input.read ( m_buffer, m_nWords+m_eventTrailerSize );
  if( !m_input.good() ){
    m_input.close();
    m_fileId++;
    if( (uint32_t)(m_fileId+1)<fileNames().size() )
      m_event=0;
    return setRunAndEventInfo(id, time, evType);
  }
  for( size_t i=0; i<m_nWords; i++ ){
    uint16_t aint;
    char buf[] = {m_buffer[i]};
    memcpy(&aint, &buf, sizeof(aint));
    rawData.push_back(aint);
    if( i<10 || i>=m_nWords-10 )
      std::cout << std::setfill('0') << std::setw(4) << std::hex << (aint) << "\t";
  }
  std::cout << std::endl;
  uint16_t evtTrailer;
  char buf[] = {m_buffer[m_nWords],m_buffer[m_nWords+1]};
  memcpy(&evtTrailer, &buf, sizeof(evtTrailer));
  std::cout << "evtTrailer = " << std::setfill('0') << std::setw(4) << std::hex << evtTrailer << std::endl;
  
  m_decodedData.clear();
  std::vector< std::array<uint16_t,1924> > decodedData=decode_raw(rawData);
  m_decodedData.insert(m_decodedData.end(),decodedData.begin(),decodedData.end());
  return true;
}

std::vector< std::array<uint16_t,1924> > HGCalTBRawDataSource::decode_raw(std::vector<uint16_t> &raw)
{
  std::vector< std::array<uint16_t, 1924> > decodedData(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA, std::array<uint16_t, 1924>());
  uint16_t x,y;
  if( !m_compressedData ){
    for(int  i = 0; i < 1924; i++){
      for (int j = 0; j < 16; j++){
	x = raw[i*16 + j];
	x = x&0xf;
	for (int sk = 0; sk < HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; sk++)
	  decodedData[sk][i] = decodedData[sk][i] | (uint16_t) (((x >> (3-sk) ) & 1) << (15 - j));
      }
    }
  }
  else{
    for(int  i = 0; i < 1924; i++){
      for (int j = 0; j < 8; j++){
	x = raw[i*8 + j];
	y = (x>>4)&0xf;
	x = x&0xf;
	for (int sk = 0; sk < HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; sk++){
	  decodedData[sk][i] = decodedData[sk][i] | (uint16_t) (((x >> (3-sk) ) & 1) << (14 - j*2));
	  decodedData[sk][i] = decodedData[sk][i] | (uint16_t) (((y >> (3-sk) ) & 1) << (15 - j*2));
	}
      }
    }
  }
  return decodedData;
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

    skirocs->push_back(skiroc);
    std::cout << std::dec << skiroc << std::endl;
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
  desc.addUntracked<unsigned int>("NumberOfBytesPerReadOut");
  desc.addUntracked<unsigned int> ("NumberOfBytesForTheHeader");
  desc.addUntracked<unsigned int> ("NumberOfBytesForTheEventTrailers");
  desc.addUntracked<unsigned int> ("NSkipEvents");
  desc.addUntracked<bool> ("CompressedData");
  descriptions.add("source", desc);
}

void HGCalTBRawDataSource::endJob()
{
}

#include "FWCore/Framework/interface/InputSourceMacros.h"
DEFINE_FWK_INPUT_SOURCE(HGCalTBRawDataSource);
