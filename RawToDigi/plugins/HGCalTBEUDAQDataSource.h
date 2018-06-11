#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include <iostream>
#include <string>
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include <TFile.h>
#include <TTree.h>
/**
 * \class HGCalTBEUDAQDataSource HGCal/RawToDigi/plugins/HGCalTBEUDAQDataSource.h
 *
 * \brief read eudaq unpacked root files and create HGCalTBSkiroc2CMS collection
 *
 */

class HGCalTBEUDAQDataSource : public edm::ProducerSourceFromFiles
{

 public:
  explicit HGCalTBEUDAQDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc);

  virtual ~HGCalTBEUDAQDataSource(){;}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);
  virtual void produce(edm::Event & e);
  virtual void endJob() override;

  bool readRaw(void);
 private:
  TFile* m_file;
  TTree* m_tree;
  
  std::string m_filePath;
  std::string m_fileName;
  std::string m_electronicMap;
  std::string m_outputCollectionName;
  unsigned int m_nSkipEvents;
  unsigned int m_fileId;
  unsigned int m_event, m_run, m_trigger;
  std::vector<uint16_t> m_rawData;
  std::vector<uint16_t> m_skiroc2CMSData;
  uint64_t m_timeStamp;

  double m_beamEnergy;
  int m_beamParticlePDGID;
  std::string m_runType;
  RUNTYPES runType;
  int m_setupConfiguration;
  

  struct {
    HGCalElectronicsMap emap_;
  } essource_;
  
  //roto ttree branches:
  int eventId;
  int triggerId;
  long long timeStamp;
  std::vector<uint16_t> *sk2cmsData;

};
