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
    float TS_MCP1;
    float TS_MCP2;
    float TS_MCP1_to_last_falling_Edge;
    float TS_MCP2_to_last_falling_Edge;

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
    tree_->Branch("TS_MCP1", &TS_MCP1);
    tree_->Branch("TS_MCP2", &TS_MCP2);
    tree_->Branch("TS_MCP1_to_last_falling_Edge", &TS_MCP1_to_last_falling_Edge);
    tree_->Branch("TS_MCP2_to_last_falling_Edge", &TS_MCP2_to_last_falling_Edge);

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

    if (rd->intUserRecords.has("valid_TS_MCP1")) valid_TS_MCP1 = rd->intUserRecords.get("valid_TS_MCP1"); else valid_TS_MCP1 = -999;
    if (rd->intUserRecords.has("valid_TS_MCP2")) valid_TS_MCP2 = rd->intUserRecords.get("valid_TS_MCP2"); else valid_TS_MCP2 = -999;
    if (rd->doubleUserRecords.has("TS_MCP1")) TS_MCP1 = rd->doubleUserRecords.get("TS_MCP1"); else TS_MCP1 = -999;
    if (rd->doubleUserRecords.has("TS_MCP2")) TS_MCP2 = rd->doubleUserRecords.get("TS_MCP2"); else TS_MCP2 = -999;
    if (rd->doubleUserRecords.has("TS_MCP1_to_last_falling_Edge")) TS_MCP1_to_last_falling_Edge = rd->doubleUserRecords.get("TS_MCP1_to_last_falling_Edge"); else TS_MCP1_to_last_falling_Edge = -999;
    if (rd->doubleUserRecords.has("TS_MCP2_to_last_falling_Edge")) TS_MCP2_to_last_falling_Edge = rd->doubleUserRecords.get("TS_MCP2_to_last_falling_Edge"); else TS_MCP2_to_last_falling_Edge = -999;

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
