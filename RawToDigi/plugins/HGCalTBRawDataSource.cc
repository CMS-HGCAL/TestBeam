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

//#define DEBUG

HGCalTBRawDataSource::HGCalTBRawDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :
  edm::ProducerSourceFromFiles(pset, desc, true),
  m_electronicMap(pset.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt")),
  m_outputCollectionName(pset.getUntrackedParameter<std::string>("OutputCollectionName","SKIROC2CMSDATA")),
  m_nWords(pset.getUntrackedParameter<unsigned int> ("NumberOf32BitsWordsPerReadOut",30788)),
  m_headerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheHeader",8)),
  m_trailerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheTrailer",4)),
  m_eventTrailerSize(pset.getUntrackedParameter<unsigned int> ("NumberOfBytesForTheEventTrailers",4)),
  m_nSkipEvents(pset.getUntrackedParameter<unsigned int> ("NSkipEvents",0)),
  m_dataFormats(pset.getUntrackedParameter<unsigned int > ("DataFormats",0)),
  m_readTimeStamps(pset.getUntrackedParameter<bool> ("ReadTimeStamps",false)),
  m_beamEnergy(pset.getUntrackedParameter<unsigned int> ("beamEnergy", 250)),
  m_beamParticlePDGID(pset.getUntrackedParameter<double> ("beamParticlePDGID", 211)),
  m_setupConfiguration(pset.getUntrackedParameter<unsigned int> ("setupConfiguration", 1))
{
  produces<HGCalTBSkiroc2CMSCollection>(m_outputCollectionName);
  produces<RunData>("RunData");
  
  m_event = 0;
  m_fileId=0;

  m_buffer=new char[m_nWords*4];//4 bytes per 32 bits
  m_header=new char[m_headerSize];
  m_trailer=new char[m_trailerSize];
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }

  m_timingFiles=pset.getParameter< std::vector<std::string> >("timingFiles");
  
  m_triggertime_prev = 0;

  std::cout << pset << std::endl;

}

bool HGCalTBRawDataSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{
  problemDuringReadout = false;

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
    if( m_readTimeStamps ) fillEventTimingInformations();
    if( m_dataFormats ){
      m_input.close();
      m_input.open(m_fileName.c_str(), std::ios::in|std::ios::binary);
    }
  }

  uint64_t nBytesToSkip=m_headerSize+m_nOrmBoards*m_nSkipEvents*(m_nWords*4+m_eventTrailerSize);
  uint32_t triggerNumber=0;
  m_decodedData.clear();
  for( uint32_t iorm=0; iorm<m_nOrmBoards; iorm++ ){
    std::vector<uint32_t> rawData;
    m_input.seekg( (std::streamoff)nBytesToSkip+(std::streamoff)m_nOrmBoards*m_event*(m_nWords*4+m_eventTrailerSize)+(std::streamoff)iorm*(m_nWords*4+m_eventTrailerSize), std::ios::beg );
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
      char buf[] = {m_buffer[i*4],m_buffer[i*4+1],m_buffer[i*4+2],m_buffer[i*4+3]};
      memcpy(&aint, &buf, sizeof(aint));
      rawData.push_back(aint);
    }
    uint32_t evtTrailer;
    char buf[] = {m_buffer[m_nWords*4],m_buffer[m_nWords*4+1],m_buffer[m_nWords*4+2],m_buffer[m_nWords*4+3]};
    memcpy(&evtTrailer, &buf, sizeof(evtTrailer));
    uint32_t ormId=evtTrailer&0xff;
    m_trigger = (evtTrailer>>0x8);
    if( ormId != iorm ) {
      problemDuringReadout = true;
      std::cout << "Problem in event trailer : wrong ORM id -> evtTrailer&0xff = " << std::dec << ormId << "\t iorm = " << iorm << std::endl;
    }
    triggerNumber=evtTrailer>>0x8;

    buf[0] = m_buffer[m_nWords*4+4];
    buf[1] = m_buffer[m_nWords*4+5];
    buf[2] = m_buffer[m_nWords*4+6];
    buf[3] = m_buffer[m_nWords*4+7];
    memcpy(&m_triggertime, &buf, sizeof(m_triggertime));
  
    std::vector< std::array<uint16_t,1924> > decodedData=decode_raw_32bit(rawData);
    m_decodedData.insert(m_decodedData.end(),decodedData.begin(),decodedData.end());
  }
  if( m_eventTimingInfoMap.find(triggerNumber)==m_eventTimingInfoMap.end() && m_readTimeStamps ){
    std::cout << "problem : cannot find trigger " << triggerNumber << " in timing information map -> exit " << std::endl;
    exit(1);
  }
  m_eventTimingInfo=m_eventTimingInfoMap[triggerNumber];

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

void HGCalTBRawDataSource::fillEventTimingInformations()
{
  switch ( m_dataFormats ){
  case 0 : readTimeStampFromTXT(); break;
  case 1  :readTimeStampFromRAW(); break;
  }
  uint64_t prevTime[m_nOrmBoards];
  uint64_t timeDiff=0;
  for(unsigned int ii=0; ii<m_nOrmBoards; ii++) prevTime[ii]=0;
  for( std::map<uint32_t,EventTimingInformation>::iterator it=m_eventTimingInfoMap.begin(); it!=m_eventTimingInfoMap.end(); ++it ){
    if( prevTime[0]==0 )
      for(unsigned int ii=0; ii<m_nOrmBoards; ii++)
	prevTime[ii]=it->second.triggerTimeStamp(ii);
    else{
      timeDiff=it->second.triggerTimeStamp(0)-prevTime[0];
      prevTime[0]=it->second.triggerTimeStamp(0);
      for(unsigned int ii=1; ii<m_nOrmBoards; ii++){
	if( it->second.triggerTimeStamp(ii)-prevTime[ii]!=timeDiff ){
	  std::cout << "Problem of timing sync in trigger " << it->first << " for rdout board " << ii
		    << " with time stamp = " << it->second.triggerTimeStamp(ii)
		    << " time diff : " << it->second.triggerTimeStamp(ii)-prevTime[ii] << " != " << timeDiff << std::endl;
        problemDuringReadout= false;
	}
	prevTime[ii]=it->second.triggerTimeStamp(ii);
      }
      //std::cout << "Trigger " << it->first << "\t time stamp = " << prevTime[0] << "\t time diff = " << timeDiff << std::endl;
    }
  }
}

void HGCalTBRawDataSource::readTimeStampFromTXT()
{
  std::ifstream input;
  std::string tc0,tc1,ts0,ts1;
  if( m_nOrmBoards != m_timingFiles.size() ) {std::cout << "oops " << m_nOrmBoards << " " << m_timingFiles.size() << std::endl; return;}
  for(uint16_t ii=0; ii<m_nOrmBoards; ii++){
    input.open(m_timingFiles[ii].c_str() , std::ifstream::in);
    if( !input.is_open() )
      std::cout << "Can't open " << m_timingFiles[ii] << "; no such file or directory"<< std::endl;
    else
      std::cout << m_timingFiles[ii] << std::endl;
    while(1){
      input >> tc0 >> tc1 >> ts0 >> ts1;
      if( !input.good() ) break; 
      char buffer[256];strcpy(buffer, tc0.c_str());
      uint32_t trig=strtoul(buffer,NULL,0);
      strcpy(buffer, ts1.c_str());
      uint64_t time=strtoul(buffer,NULL,0);
      if( m_eventTimingInfoMap.find(trig)==m_eventTimingInfoMap.end() ){
	EventTimingInformation ETI(trig);
	ETI.addTriggerTimeStamp(time);
	std::pair<uint32_t,EventTimingInformation> p(trig,ETI);
	m_eventTimingInfoMap.insert(p);
      }
      else
	m_eventTimingInfoMap[trig].addTriggerTimeStamp(time);
    }
    input.close();
  }
}

void HGCalTBRawDataSource::readTimeStampFromRAW()
{
  uint64_t nBytesToSkip=m_headerSize+m_nOrmBoards*m_nSkipEvents*(m_nWords*4+m_eventTrailerSize);
  uint32_t triggerNumber=0;
  int nevent=0;
  bool letsbreak=false;
  while(1){
    for( uint32_t iorm=0; iorm<m_nOrmBoards; iorm++ ){
      m_input.seekg( (std::streamoff)nBytesToSkip+(std::streamoff)m_nOrmBoards*nevent*(m_nWords*4+m_eventTrailerSize)+(std::streamoff)iorm*(m_nWords*4+m_eventTrailerSize), std::ios::beg );
      m_input.read ( m_buffer, m_nWords*4+m_eventTrailerSize );
      if( !m_input.good() ) {letsbreak=true;break;}
      uint32_t evtTrailer[3];
      char buf[] = {m_buffer[m_nWords*4],m_buffer[m_nWords*4+1],m_buffer[m_nWords*4+2],m_buffer[m_nWords*4+3]};
      memcpy(&evtTrailer[0], &buf, sizeof(evtTrailer[0]));
      char ts0[] = {m_buffer[m_nWords*4+4],m_buffer[m_nWords*4+5],m_buffer[m_nWords*4+6],m_buffer[m_nWords*4+7]};
      memcpy(&evtTrailer[1], &ts0, sizeof(evtTrailer[1]));
      char ts1[] = {m_buffer[m_nWords*4+8],m_buffer[m_nWords*4+9],m_buffer[m_nWords*4+10],m_buffer[m_nWords*4+11]};//event trailer size must be at least 12 when dataformats = 1
      memcpy(&evtTrailer[2], &ts1, sizeof(evtTrailer[2]));
    
      triggerNumber=evtTrailer[0]>>0x8;    
      uint64_t timeStamp=((uint64_t)evtTrailer[2]<<32)|(uint64_t)evtTrailer[1];
      if( m_eventTimingInfoMap.find(triggerNumber)==m_eventTimingInfoMap.end() ){
	EventTimingInformation ETI(triggerNumber);
	ETI.addTriggerTimeStamp(timeStamp);
	std::pair<uint32_t,EventTimingInformation> p(triggerNumber,ETI);
	m_eventTimingInfoMap.insert(p);
      }
      else
	m_eventTimingInfoMap[triggerNumber].addTriggerTimeStamp(timeStamp);
    }
    if( letsbreak ) break;
    nevent++;
  }
}

void HGCalTBRawDataSource::produce(edm::Event & e)
{
  std::auto_ptr<HGCalTBSkiroc2CMSCollection> skirocs(new HGCalTBSkiroc2CMSCollection);
  std::auto_ptr<HGCalTBSkiroc2CMSCollection> emptycol(new HGCalTBSkiroc2CMSCollection);
  //bool correctEvent=true;
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
    HGCalTBSkiroc2CMS skiroc;
    if( m_readTimeStamps )
      skiroc=HGCalTBSkiroc2CMS( vdata,detids,
				m_eventTimingInfo.triggerTimeStamp(0),
				m_eventTimingInfo.triggerCounter());
    else skiroc=HGCalTBSkiroc2CMS( vdata,detids );
    if(skiroc.check())
      skirocs->push_back(skiroc);
    else{
      std::cout<<"Skiroc info: "<<iski<<"  "<<iski/4<<std::endl;
      //correctEvent=false;
      //break;
    }
  }
  //if( !correctEvent )
    //skirocs->swap(*emptycol);
  e.put(skirocs, m_outputCollectionName);


  m_event++;
  

  //set the RunData
  std::auto_ptr<RunData> rd(new RunData);


  rd->configuration = m_setupConfiguration;
  rd->energy = m_beamEnergy;
  rd->runType = "HGCal_TB";
  rd->pdgID = m_beamParticlePDGID;
  rd->run = m_run;
  rd->event = m_event;
  rd->booleanUserRecords.add("hasDanger", problemDuringReadout);

  #ifdef DEBUG
    std::cout<<rd->run<<"  "<<rd->event<<"  "<<rd->energy<<"  "<<rd->configuration<<"  "<<rd->runType<<"  "<<rd->booleanUserRecords.get("hasDanger")<<std::endl;
  #endif

  e.put(std::move(rd), "RunData");
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
  desc.addUntracked<unsigned int> ("NSkipEvents");
  desc.addUntracked<bool> ("ReadTimeStamps");
  desc.addUntracked<unsigned int> ("DataFormats");
  desc.add<std::vector<std::string> >("timingFiles");
  desc.addUntracked<unsigned int> ("beamEnergy");
  desc.addUntracked<double> ("beamParticlePDGID");
  desc.addUntracked<unsigned int> ("setupConfiguration");

  descriptions.add("source", desc);
}

void HGCalTBRawDataSource::endJob()
{
  // std::cout << "mean reading = " << m_meanReadingTime/m_event
  // 	    << "  +-  " << std::sqrt(m_rmsReadingTime/m_event-m_meanReadingTime/m_event*m_meanReadingTime/m_event)
  // 	    << std::endl;
}

#include "FWCore/Framework/interface/InputSourceMacros.h"
DEFINE_FWK_INPUT_SOURCE(HGCalTBRawDataSource);
