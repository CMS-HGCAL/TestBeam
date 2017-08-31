#include <iostream>
#include <fstream>
#include <sstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMSCollection.h"

#include <iomanip>
#include <set>
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

class HGCalTBTimingFileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit HGCalTBTimingFileWriter(const edm::ParameterSet&);
  ~HGCalTBTimingFileWriter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;
  edm::EDGetTokenT<HGCalTBSkiroc2CMSCollection> m_HGCalTBSkiroc2CMSCollection;
  std::string m_TimingFilePath;
  std::ofstream timingFile;
  uint64_t m_prev_time;
  int m_event;

};

HGCalTBTimingFileWriter::HGCalTBTimingFileWriter(const edm::ParameterSet& iConfig) :
  m_TimingFilePath(iConfig.getUntrackedParameter<std::string>("TimingFilePath",""))
{
  m_HGCalTBSkiroc2CMSCollection = consumes<HGCalTBSkiroc2CMSCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));


}


HGCalTBTimingFileWriter::~HGCalTBTimingFileWriter() {

}

void HGCalTBTimingFileWriter::beginJob() {
    timingFile.open(m_TimingFilePath);
    if (timingFile.is_open()) {
      timingFile<<"TrigNumber"<<"   "<<"TrigCount"<<"   "<<"TimeStamp"<<"   "<<"TimeDiff"<<std::endl;
    }
    m_event=0;
    m_prev_time=0;
}


void HGCalTBTimingFileWriter::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  m_event++;

  edm::Handle<HGCalTBSkiroc2CMSCollection> skirocs;
  event.getByToken(m_HGCalTBSkiroc2CMSCollection, skirocs);

  uint64_t triggerTime = skirocs->at(0).triggerTimeStamp();
  int triggerCount = skirocs->at(0).triggerCounter();

    if (timingFile.is_open()) {
      timingFile<<m_event<<"   "<<triggerCount<<"   "<<triggerTime<<"   "<<triggerTime-m_prev_time<<std::endl;
      std::cout<<m_event<<"   "<<triggerCount<<"   "<<triggerTime<<"   "<<triggerTime-m_prev_time<<std::endl;
    }

  m_prev_time = triggerTime;

}

void HGCalTBTimingFileWriter::endJob()
{
}

void HGCalTBTimingFileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HGCalTBTimingFileWriter);
