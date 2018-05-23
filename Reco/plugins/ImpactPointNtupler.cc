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
#include "HGCal/DataFormats/interface/HGCalTBDATURATelescopeData.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h" //for the runData type definition

#include <iomanip>
#include <set>

class ImpactPointNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit ImpactPointNtupler(const edm::ParameterSet&);
    ~ImpactPointNtupler();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void endJob() override;
    void analyze(const edm::Event& , const edm::EventSetup&) override;

    

    // ---------- member data ---------------------------
    edm::EDGetTokenT<std::vector<HGCalTBDATURATelescopeData> > DATURATrackToken;    
    edm::EDGetTokenT<RunData> RunDataToken;

    // Output tree
    TTree* tree_;

    void clearVariables(); // function to clear tree variables/vectors

    // Variables for branches

    // event info
    unsigned int ev_run_;
    unsigned int ev_event_;

    int nTrackCounter;
    // impact points
    std::map<int, std::vector<float> >impactX;
    std::map<int, std::vector<float> >impactY;
    std::map<int, std::vector<float> >impactX_associatedChi2;
    std::map<int, std::vector<float> >impactY_associatedChi2;
    std::vector<float> kinkAngleX_DUT1;
    std::vector<float> kinkAngleY_DUT1;


    int m_nLayers;
};

void ImpactPointNtupler::clearVariables(){
    // event info
    ev_run_ = 0;
    ev_event_ = 0;


    for (int layer=1; layer<=m_nLayers; layer++) {
        impactX[layer].clear();
        impactY[layer].clear();
        impactX_associatedChi2[layer].clear();
        impactY_associatedChi2[layer].clear();        
    }
    kinkAngleX_DUT1.clear();
    kinkAngleY_DUT1.clear();    
};

ImpactPointNtupler::ImpactPointNtupler(const edm::ParameterSet& iConfig) 
{

    DATURATrackToken= consumes<std::vector<HGCalTBDATURATelescopeData> >(iConfig.getParameter<edm::InputTag>("DATURATelescopeData"));
    RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
    m_nLayers= iConfig.getUntrackedParameter<int>("nLayers", 10);


    for (int layer=1; layer<=m_nLayers; layer++) {
        impactX[layer] = std::vector<float>(0);
        impactY[layer] = std::vector<float>(0);
        impactX_associatedChi2[layer] = std::vector<float>(0);
        impactY_associatedChi2[layer] = std::vector<float>(0);
    }

    usesResource("TFileService");
    edm::Service<TFileService> fs;

    // Define tree and branches
    tree_ = fs->make<TTree>("impactPoints", "impactPoints");

    // event info
    tree_->Branch("event", &ev_event_);
    tree_->Branch("run", &ev_run_);
    tree_->Branch("ntracks", &nTrackCounter);
    for (int layer=1; layer<=m_nLayers; layer++) {
        tree_->Branch(("impactX_HGCal_layer_"+std::to_string(layer)).c_str(), &impactX[layer]);
        tree_->Branch(("impactY_HGCal_layer_"+std::to_string(layer)).c_str(), &impactY[layer]);
        tree_->Branch(("impactX_associatedChi2_HGCal_layer_"+std::to_string(layer)).c_str(), &impactX_associatedChi2[layer]);
        tree_->Branch(("impactY_associatedChi2_HGCal_layer_"+std::to_string(layer)).c_str(), &impactY_associatedChi2[layer]);
    }  
    tree_->Branch("kinkAngleX_DUT1", &kinkAngleX_DUT1);  
    tree_->Branch("kinkAngleY_DUT1", &kinkAngleY_DUT1);  
}


ImpactPointNtupler::~ImpactPointNtupler()
{

}

void ImpactPointNtupler::beginJob()
{
}

void ImpactPointNtupler::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    clearVariables();

    edm::Handle<RunData> rd;
    event.getByToken(RunDataToken, rd);    
    ev_run_ = rd->run;
    ev_event_ = rd->event;


    edm::Handle<std::vector<HGCalTBDATURATelescopeData> > daturatracks; 
    event.getByToken(DATURATrackToken, daturatracks);
    nTrackCounter=0;
    
    for(auto daturatrack : *daturatracks) {
        nTrackCounter++;
        for(int layer=1; layer<=m_nLayers; layer++) {        
            if (daturatrack.Extrapolation_XY_Chi2(layer).first > 100.) continue;
            if (daturatrack.Extrapolation_XY_Chi2(layer).second > 100.) continue;
            daturatrack.Extrapolation_XY(layer).first;

            impactX[layer].push_back(daturatrack.Extrapolation_XY(layer).first);
            impactY[layer].push_back(daturatrack.Extrapolation_XY(layer).second);
            impactX_associatedChi2[layer].push_back(daturatrack.Extrapolation_XY_Chi2(layer).first);
            impactY_associatedChi2[layer].push_back(daturatrack.Extrapolation_XY_Chi2(layer).second);
        }

        if(daturatrack.floatUserRecords.has("kinkAngleX_DUT1")) kinkAngleX_DUT1.push_back(daturatrack.floatUserRecords.get("kinkAngleX_DUT1"));
        if(daturatrack.floatUserRecords.has("kinkAngleY_DUT1")) kinkAngleY_DUT1.push_back(daturatrack.floatUserRecords.get("kinkAngleY_DUT1"));
    }
    
    tree_->Fill();
}


void ImpactPointNtupler::endJob()
{
}

void ImpactPointNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ImpactPointNtupler);
