#include "HGCal/RawToDigi/plugins/HGCalTBEUDAQDataSource.h"
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

HGCalTBEUDAQDataSource::HGCalTBEUDAQDataSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :
  edm::ProducerSourceFromFiles(pset, desc, true),
  m_electronicMap(pset.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers.txt")),
  m_outputCollectionName(pset.getUntrackedParameter<std::string>("OutputCollectionName","SKIROC2CMSDATA")),
  m_nSkipEvents(pset.getUntrackedParameter<unsigned int> ("NSkipEvents",0)),
  m_beamEnergy(pset.getUntrackedParameter<double> ("beamEnergy", 0)),
  m_beamParticlePDGID(pset.getUntrackedParameter<int> ("beamParticlePDGID", 0)),
  m_runType(pset.getUntrackedParameter<std::string> ("runType", "Beam")),
  m_setupConfiguration(pset.getUntrackedParameter<unsigned int> ("setupConfiguration", 1))
{
  produces<HGCalTBSkiroc2CMSCollection>(m_outputCollectionName);
  produces<RunData>("RunData");
  
  m_event = 0;
  m_fileId=0;
  m_skiroc2CMSData.resize(1924);

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(m_electronicMap);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }

  std::cout << pset << std::endl;

  if (m_runType=="Pedestal") {
    runType = HGCAL_TB_PEDESTAL;
  } else if (m_runType=="Beam") {
    runType = HGCAL_TB_BEAM;
  } else if (m_runType=="Simulation") {
    runType = HGCAL_TB_PEDESTAL;
  } else {
    runType = HGCAL_TB_BEAM;
  }

}

bool HGCalTBEUDAQDataSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{
  if( m_fileId == fileNames().size() ) return false;
  if (fileNames()[m_fileId] != "file:DUMMY")
    m_fileName = fileNames()[m_fileId];
  if (m_fileName.find("file:") == 0) m_fileName = m_fileName.substr(5);
  if( m_event==0 ){
    m_file = new TFile(m_fileName.c_str(),"READ");
    if( m_file->IsOpen() ){
      m_file->Print();
    }

    else{
      std::cout << "can not open file " << m_fileName << std::endl;
      return false;
    }
    m_tree = (TTree*)m_file->Get("hgcalraw");
    if (!m_tree){
      std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
      return 0;
    }
    m_tree->SetBranchAddress("eventId", &eventId);
    m_tree->SetBranchAddress("triggerId", &triggerId);
    m_tree->SetBranchAddress("timeStamp", &timeStamp);
    m_tree->SetBranchAddress("sk2cmsData", &sk2cmsData);
    sk2cmsData=0;
  }

  if( m_event>=(uint32_t)m_tree->GetEntries() ){
    m_fileId++;
    if( (uint32_t)(m_fileId+1)<fileNames().size() )
      m_event=0;
    return setRunAndEventInfo(id, time, evType);
  }

  m_tree->GetEntry(m_event);
  if(m_event!=(uint32_t)eventId){
    std::cout << "problem in event ID : m_event = " << m_event << "; eventId = " << eventId << std::endl;
    return false;
  }
  if(m_event==0)
    m_rawData.resize( sk2cmsData->size() );
  std::copy(sk2cmsData->begin(),sk2cmsData->end(),m_rawData.begin());
  m_timeStamp=timeStamp;
  m_trigger=triggerId;
  m_event++;
  return true;
}

void HGCalTBEUDAQDataSource::produce(edm::Event & e)
{
  std::unique_ptr<HGCalTBSkiroc2CMSCollection> skirocs(new HGCalTBSkiroc2CMSCollection);

  uint32_t vectorIndex=0;
  unsigned int skirocId=0;
  while(1){
    if(vectorIndex>=m_rawData.size()) break;
    std::vector<HGCalTBDetId> detids;
    int skiId=HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA*(skirocId/HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA)+(HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA-skirocId)%HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA+1;
    for (size_t ichan = 0; ichan < HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++) {
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
    std::copy(m_rawData.begin()+vectorIndex,m_rawData.begin()+vectorIndex+1924,m_skiroc2CMSData.begin());
    HGCalTBSkiroc2CMS skiroc;
    skiroc=HGCalTBSkiroc2CMS( m_skiroc2CMSData,detids,
			      m_timeStamp,
			      m_trigger);
    skirocs->push_back(skiroc);
    vectorIndex+=1924;
    skirocId++;
  }
  e.put(std::move(skirocs), m_outputCollectionName);


    
  //set the RunData
  std::unique_ptr<RunData> rd(new RunData);

  rd->configuration = m_setupConfiguration;
  rd->energy = m_beamEnergy;
  rd->runType = runType;
  rd->pdgID = m_beamParticlePDGID;
  rd->run = m_run;
  rd->event = m_event;
  //rd->booleanUserRecords.add("hasDanger", problemDuringReadout);

  #ifdef DEBUG
    std::cout<<rd->run<<"  "<<rd->event<<"  "<<rd->energy<<"  "<<rd->configuration<<"  "<<rd->runType<<"  "<<rd->booleanUserRecords.get("hasDanger")<<std::endl;
  #endif

  e.put(std::move(rd), "RunData");
}

void HGCalTBEUDAQDataSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setComment("TEST");
  desc.addUntracked<int>("run", 101);
  desc.addUntracked<std::vector<std::string> >("fileNames");
  desc.addUntracked<std::string>("ElectronicMap");
  desc.addUntracked<std::string>("OutputCollectionName");
  desc.addUntracked<unsigned int> ("NSkipEvents");
  desc.addUntracked<double> ("beamEnergy");
  desc.addUntracked<int> ("beamParticlePDGID");
  desc.addUntracked<std::string>("runType");
  desc.addUntracked<unsigned int> ("setupConfiguration");

  descriptions.add("source", desc);
}

void HGCalTBEUDAQDataSource::endJob()
{
  // std::cout << "mean reading = " << m_meanReadingTime/m_event
  // 	    << "  +-  " << std::sqrt(m_rmsReadingTime/m_event-m_meanReadingTime/m_event*m_meanReadingTime/m_event)
  // 	    << std::endl;
}

#include "FWCore/Framework/interface/InputSourceMacros.h"
DEFINE_FWK_INPUT_SOURCE(HGCalTBEUDAQDataSource);
