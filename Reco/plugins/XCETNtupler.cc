#include <iostream>
#include "TTree.h"
#include <fstream>
#include <sstream>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h" //for the runData type definition

#include <iomanip>
#include <set>



class XCETNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit XCETNtupler(const edm::ParameterSet&);
    ~XCETNtupler();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void endJob() override;
    void analyze(const edm::Event& , const edm::EventSetup&) override;



    // ---------- member data ---------------------------
    edm::EDGetTokenT<RunData> RunDataToken;


    // Output tree
    TTree* tree_;

    void clearVariables(); // function to clear tree variables/vectors

    // Variables for branches

    // event info
    unsigned int ev_run_;
    unsigned int ev_event_;

    short XCET_021507_signal;
    short XCET_021523_signal;
    short scintillator_coincidences;
    short scintillator_vetos;

};

void XCETNtupler::clearVariables() {
    // event info
    ev_run_ = 0;
    ev_event_ = 0;
}

XCETNtupler::XCETNtupler(const edm::ParameterSet& iConfig)
{
    RunDataToken = consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));


    usesResource("TFileService");
    edm::Service<TFileService> fs;

    // Define tree and branches
    tree_ = fs->make<TTree>("XCET", "XCET");

    // event info
    tree_->Branch("event", &ev_event_);
    tree_->Branch("run", &ev_run_);
    tree_->Branch("XCET_021507_signal", &XCET_021507_signal);
    tree_->Branch("XCET_021523_signal", &XCET_021523_signal);
    tree_->Branch("scintillator_coincidences", &scintillator_coincidences);
    tree_->Branch("scintillator_vetos", &scintillator_vetos);

}


XCETNtupler::~XCETNtupler()
{

}

void XCETNtupler::beginJob()
{
}

void XCETNtupler::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    clearVariables();

    edm::Handle<RunData> rd;
    event.getByToken(RunDataToken, rd);
    ev_run_ = rd->run;
    ev_event_ = rd->event;

    if (rd->intUserRecords.has("XCET_021507_signal")) XCET_021507_signal = rd->intUserRecords.get("XCET_021507_signal"); else XCET_021507_signal = -999;
    if (rd->intUserRecords.has("XCET_021523_signal")) XCET_021523_signal = rd->intUserRecords.get("XCET_021523_signal"); else XCET_021523_signal = -999;
    if (rd->intUserRecords.has("scintillator_coincidences")) scintillator_coincidences = rd->intUserRecords.get("scintillator_coincidences"); else scintillator_coincidences = -999;
    if (rd->intUserRecords.has("scintillator_vetos")) scintillator_vetos = rd->intUserRecords.get("scintillator_vetos"); else scintillator_vetos = -999;

    tree_->Fill();
}


void XCETNtupler::endJob()
{
}

void XCETNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(XCETNtupler);
