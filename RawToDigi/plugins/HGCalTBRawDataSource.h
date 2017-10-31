#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "stdlib.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
/**
 * \class HGCalTBRawDataSource HGCal/RawToDigi/plugins/HGCalTBRawDataSource.h
 *
 * \brief convert data from raw file to HGCalTBSkiroc2CMS data
 *
 */

class EventTimingInformation{
 public:
  EventTimingInformation(){;}
  EventTimingInformation(const uint32_t trigger){m_triggerCounter=trigger;}
  uint32_t triggerCounter() const {return m_triggerCounter;}
  uint64_t triggerTimeStamp(int ormId) const {return m_triggerTimeStamp.at(ormId);}
  void addTriggerTimeStamp( uint64_t ts ){ m_triggerTimeStamp.push_back(ts); }
 private:
  uint32_t m_triggerCounter;
  std::vector<uint64_t> m_triggerTimeStamp;
};

class HGCalTBRawDataSource : public edm::ProducerSourceFromFiles
{

 public:
  explicit HGCalTBRawDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);

  virtual ~HGCalTBRawDataSource(){;}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);
  virtual void produce(edm::Event & e);
  virtual void endJob() override;

  bool readRaw(void);
  std::vector< std::array<uint16_t,1924> > decode_raw_32bit(std::vector<uint32_t> &raw);
  void fillEventTimingInformations();
 private:
  std::string m_filePath;
  std::string m_fileName;
  std::string m_electronicMap;
  std::string m_outputCollectionName;
  std::vector<std::string> m_timingFiles;
  //std::vector<FileInfo>::iterator fileIterator;
  unsigned int m_nHexaboards;
  unsigned int m_nWords;
  unsigned int m_headerSize;
  unsigned int m_trailerSize;
  unsigned int m_eventTrailerSize;
  unsigned int m_nOrmBoards;
  unsigned int m_nSkipEvents;
  unsigned int m_dataFormats; //0 when time stamps were save in separated .txt files, 1 since the timestamp are saved in the raw data
  bool m_readTimeStamps;

  int m_beamEnergy;
  std::string m_beamParticlePDGID;
  int m_setupConfiguration;
  
  unsigned int timeStartRun;
  unsigned int timeStopRun;
  
  unsigned int m_time;
  unsigned int m_event, m_run, m_spill, m_trigger, m_triggertime, m_triggertime_prev;
  int m_eventSize,m_nevents;
  unsigned int m_fileId;
  
  bool problemDuringReadout;
  
  uint32_t m_skiMask;
  char* m_buffer;
  char* m_header;
  char* m_trailer;
  std::ifstream m_input;
  std::vector< std::array<uint16_t,1924> > m_decodedData;
  std::map<uint32_t,EventTimingInformation> m_eventTimingInfoMap;
  EventTimingInformation m_eventTimingInfo;
  
  struct {
    HGCalElectronicsMap emap_;
  } essource_;

  void readTimeStampFromTXT();
  void readTimeStampFromRAW();
};
