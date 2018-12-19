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



class MCPNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit MCPNtupler(const edm::ParameterSet&);
    ~MCPNtupler();
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

    short valid_TS_MCP1;
    short valid_TS_MCP2;
    float noise_MCP1;
    float noise_MCP2;
    float TSpeak_MCP1;
    float TSpeak_MCP2;
    float amp_MCP1;
    float amp_MCP2;
    float ampFit_MCP1;
    float ampFit_MCP2;
    float TSfitPeak_MCP1;
    float TSfitPeak_MCP2;
    float TScf_MCP1;
    float TScf_MCP2;
    float charge5nsS_MCP1;
    float charge5nsS_MCP2;
    float charge5nsB_MCP1;
    float charge5nsB_MCP2;
    float TS_toClock_FE_MCP1;
    float TS_toClock_FE_MCP2;
    float meanClockFE;
    float rmsClockFE;
};

void MCPNtupler::clearVariables() {
    // event info
    ev_run_ = 0;
    ev_event_ = 0;
}

MCPNtupler::MCPNtupler(const edm::ParameterSet& iConfig)
{
    RunDataToken = consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));


    usesResource("TFileService");
    edm::Service<TFileService> fs;

    // Define tree and branches
    tree_ = fs->make<TTree>("MCP", "MCP");

    // event info
    tree_->Branch("event", &ev_event_);
    tree_->Branch("run", &ev_run_);

    tree_->Branch("valid_TS_MCP1", &valid_TS_MCP1);
    tree_->Branch("valid_TS_MCP2", &valid_TS_MCP2);
    tree_->Branch("noise_MCP1", &noise_MCP1);
    tree_->Branch("noise_MCP2", &noise_MCP2);
    tree_->Branch("TSpeak_MCP1", &TSpeak_MCP1); 
    tree_->Branch("TSpeak_MCP2", &TSpeak_MCP2); 
    tree_->Branch("amp_MCP1", &amp_MCP1);
    tree_->Branch("amp_MCP2", &amp_MCP2);
    tree_->Branch("ampFit_MCP1", &ampFit_MCP1);
    tree_->Branch("ampFit_MCP2", &ampFit_MCP2);
    tree_->Branch("TSfitPeak_MCP1", &TSfitPeak_MCP1);
    tree_->Branch("TSfitPeak_MCP2", &TSfitPeak_MCP2);
    tree_->Branch("TScf_MCP1", &TScf_MCP1);
    tree_->Branch("TScf_MCP2", &TScf_MCP2);
    tree_->Branch("charge5nsS_MCP1", &charge5nsS_MCP1);
    tree_->Branch("charge5nsS_MCP2", &charge5nsS_MCP2);
    tree_->Branch("charge5nsB_MCP1", &charge5nsB_MCP1);
    tree_->Branch("charge5nsB_MCP2", &charge5nsB_MCP2);
    tree_->Branch("TS_toClock_FE_MCP1", &TS_toClock_FE_MCP1);
    tree_->Branch("TS_toClock_FE_MCP2", &TS_toClock_FE_MCP2);
    tree_->Branch("meanClockFE", &meanClockFE);
    tree_->Branch("rmsClockFE", &rmsClockFE);
}


MCPNtupler::~MCPNtupler()
{

}

void MCPNtupler::beginJob()
{
}

void MCPNtupler::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    clearVariables();

    edm::Handle<RunData> rd;
    event.getByToken(RunDataToken, rd);
    ev_run_ = rd->run;
    ev_event_ = rd->event;


    valid_TS_MCP1 = (rd->intUserRecords.has("valid_TS_MCP1")) ? rd->intUserRecords.get("valid_TS_MCP1") : -999;
    valid_TS_MCP2 = (rd->intUserRecords.has("valid_TS_MCP2")) ? rd->intUserRecords.get("valid_TS_MCP2") : -999;
    noise_MCP1 = (rd->doubleUserRecords.has("noise_MCP1")) ? rd->doubleUserRecords.get("noise_MCP1") : -999;
    noise_MCP2 = (rd->doubleUserRecords.has("noise_MCP2")) ? rd->doubleUserRecords.get("noise_MCP2") : -999;
    TSpeak_MCP1 = (rd->doubleUserRecords.has("TSpeak_MCP1")) ? rd->doubleUserRecords.get("TSpeak_MCP1") : -999;
    TSpeak_MCP2 = (rd->doubleUserRecords.has("TSpeak_MCP2")) ? rd->doubleUserRecords.get("TSpeak_MCP2") : -999;
    amp_MCP1 = (rd->doubleUserRecords.has("amp_MCP1")) ? rd->doubleUserRecords.get("amp_MCP1") : -999;
    amp_MCP2 = (rd->doubleUserRecords.has("amp_MCP2")) ? rd->doubleUserRecords.get("amp_MCP2") : -999;
    ampFit_MCP1 = (rd->doubleUserRecords.has("ampFit_MCP1")) ? rd->doubleUserRecords.get("ampFit_MCP1") : -999;
    ampFit_MCP2 = (rd->doubleUserRecords.has("ampFit_MCP2")) ? rd->doubleUserRecords.get("ampFit_MCP2") : -999;
    TSfitPeak_MCP1 = (rd->doubleUserRecords.has("TSfitPeak_MCP1")) ? rd->doubleUserRecords.get("TSfitPeak_MCP1") : -999;
    TSfitPeak_MCP2 = (rd->doubleUserRecords.has("TSfitPeak_MCP2")) ? rd->doubleUserRecords.get("TSfitPeak_MCP2") : -999;
    TScf_MCP1 = (rd->doubleUserRecords.has("TScf_MCP1")) ? rd->doubleUserRecords.get("TScf_MCP1") : -999;
    TScf_MCP2 = (rd->doubleUserRecords.has("TScf_MCP2")) ? rd->doubleUserRecords.get("TScf_MCP2") : -999;
    charge5nsS_MCP1 = (rd->doubleUserRecords.has("charge5nsS_MCP1")) ? rd->doubleUserRecords.get("charge5nsS_MCP1") : -999;
    charge5nsS_MCP2 = (rd->doubleUserRecords.has("charge5nsS_MCP2")) ? rd->doubleUserRecords.get("charge5nsS_MCP2") : -999;
    charge5nsB_MCP1 = (rd->doubleUserRecords.has("charge5nsB_MCP1")) ? rd->doubleUserRecords.get("charge5nsB_MCP1") : -999;
    charge5nsB_MCP2 = (rd->doubleUserRecords.has("charge5nsB_MCP2")) ? rd->doubleUserRecords.get("charge5nsB_MCP2") : -999;
    TS_toClock_FE_MCP1 = (rd->doubleUserRecords.has("TS_toClock_FE_MCP1")) ? rd->doubleUserRecords.get("TS_toClock_FE_MCP1") : -999;
    TS_toClock_FE_MCP2 = (rd->doubleUserRecords.has("TS_toClock_FE_MCP2")) ? rd->doubleUserRecords.get("TS_toClock_FE_MCP2") : -999;
    meanClockFE = (rd->doubleUserRecords.has("meanClockFE")) ? rd->doubleUserRecords.get("meanClockFE") : -999;
    rmsClockFE = (rd->doubleUserRecords.has("rmsClockFE")) ? rd->doubleUserRecords.get("rmsClockFE") : -999;

    tree_->Fill();
}


void MCPNtupler::endJob()
{
}

void MCPNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MCPNtupler);
