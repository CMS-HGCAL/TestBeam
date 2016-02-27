#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include <stdio.h>

class HGCalTBTextSource : public edm::ProducerSourceFromFiles {
public:
  explicit HGCalTBTextSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset,desc, true),
													 m_file(0),
													 m_run(pset.getUntrackedParameter<int>("run",101))
																															
  {
    
    m_sourceId=pset.getUntrackedParameter<int>("fed",1000);
    produces<FEDRawDataCollection>();    
  }
  virtual ~HGCalTBTextSource() {
    if (m_file!=0) fclose(m_file); m_file=0;
  }
    
private:
  virtual bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&) {

    if (!(fileNames().size())) return false; // need a file...
    if (m_file==0) {
      m_file=fopen(fileNames()[0].c_str(),"r");
      if (m_file==0) return false; // can't open the file!
      m_event=0;
    }
    if (feof(m_file)) return false;
    
    m_event++;
    id=edm::EventID(m_run,1,m_event);
    // time is a hack
    edm::TimeValue_t present_time = presentTime();
    unsigned long time_between_events = timeBetweenEvents();
    
    time = present_time + time_between_events;
    return true;
  }
  virtual void produce(edm::Event & e);

  FILE* m_file;
  int m_event, m_run;
  int m_sourceId;
};

void HGCalTBTextSource::produce(edm::Event & e) {

}


#include "FWCore/Framework/interface/InputSourceMacros.h"
 
DEFINE_FWK_INPUT_SOURCE(HGCalTBTextSource);
