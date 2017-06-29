#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
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
 * \brief convert data from raw file to FEDRawData
 *
 */


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

 private:
  std::string m_filePath;
  std::string m_fileName;
  std::string m_electronicMap;
  std::string m_outputCollectionName;
  //std::vector<FileInfo>::iterator fileIterator;
  unsigned int m_nOrmBoards;
  unsigned int m_nHexaboards;
  unsigned int m_nWords;

  unsigned int m_time;
  unsigned int m_event, m_run, m_spill;
  int m_eventSize,m_nevents;
  unsigned int m_fileId;
  
  uint32_t m_skiMask;
  char* m_buffer;
  std::ifstream m_input;
  std::vector< std::array<uint16_t,1924> > m_decodedData;
  
  struct {
    HGCalElectronicsMap emap_;
  } essource_;
};
