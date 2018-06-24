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
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h" //for the runData type definition

#include <iomanip>
#include <set>

enum ExtrapolationType {
    DATURA = 0,
    DWC 
};

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
    edm::EDGetTokenT<HGCalTBDWCTrack> DWCTrackToken;    
    edm::EDGetTokenT<RunData> RunDataToken;

    //DWC track or DATURA track
    std::string extrapolationTypeString;
    ExtrapolationType _extrapolationType;

    // Output tree
    TTree* tree_;

    void clearVariables(); // function to clear tree variables/vectors

    // Variables for branches

    // event info
    unsigned int ev_run_;
    unsigned int ev_event_;

    int nTrackCounter;
    // impact points
    std::map<int, std::vector<float> >impactsX;
    std::map<int, std::vector<float> >impactsY;
    std::map<int, std::vector<float> >impactsX_associatedChi2;
    std::map<int, std::vector<float> >impactsY_associatedChi2;

    std::map<int, float>impactX;
    std::map<int, float>impactY;
    std::map<int, float>impactX_associatedChi2;
    std::map<int, float>impactY_associatedChi2;
    
    //specific to 6 plane DATURA
    std::vector<float> kinkAngleX_DUT1;
    std::vector<float> kinkAngleY_DUT1;

    //specific to DWC in H2
    int dwcReferenceType;
    double m_x, m_y; //slopes of straight line tracks in x/y
    double b_x, b_y; //offsets of straight line tracks in x/y

    int m_nLayers;
};

void ImpactPointNtupler::clearVariables(){
    // event info
    ev_run_ = 0;
    ev_event_ = 0;


    for (int layer=0; layer<=m_nLayers; layer++) {
        impactsX[layer].clear();
        impactsY[layer].clear();
        impactsX_associatedChi2[layer].clear();
        impactsY_associatedChi2[layer].clear();                
    }
    kinkAngleX_DUT1.clear();
    kinkAngleY_DUT1.clear();    
};

ImpactPointNtupler::ImpactPointNtupler(const edm::ParameterSet& iConfig) 
{
    extrapolationTypeString= iConfig.getUntrackedParameter<std::string>("extrapolationDevice", "DATURA");
    if (extrapolationTypeString=="DATURA") _extrapolationType = DATURA;
    else if (extrapolationTypeString=="DWC") _extrapolationType = DWC;
    else _extrapolationType = DATURA;

    DATURATrackToken= consumes<std::vector<HGCalTBDATURATelescopeData> >(iConfig.getParameter<edm::InputTag>("DATURATelescopeData"));
    DWCTrackToken= consumes<HGCalTBDWCTrack >(iConfig.getParameter<edm::InputTag>("DWCTrackToken"));
    RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
    m_nLayers= iConfig.getUntrackedParameter<int>("nLayers", 10);


    for (int layer=1; layer<=m_nLayers; layer++) {
        impactsX[layer] = std::vector<float>(0);
        impactsY[layer] = std::vector<float>(0);
        impactsX_associatedChi2[layer] = std::vector<float>(0);
        impactsY_associatedChi2[layer] = std::vector<float>(0);
        
        impactX[layer] = -999;
        impactY[layer] = -999;
        impactX_associatedChi2[layer] = -999;
        impactY_associatedChi2[layer] = -999;
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
        if (_extrapolationType == DATURA) {
            tree_->Branch(("impactX_HGCal_layer_"+std::to_string(layer)).c_str(), &impactsX[layer]);
            tree_->Branch(("impactY_HGCal_layer_"+std::to_string(layer)).c_str(), &impactsY[layer]);
            tree_->Branch(("impactX_associatedChi2_HGCal_layer_"+std::to_string(layer)).c_str(), &impactsX_associatedChi2[layer]);
            tree_->Branch(("impactY_associatedChi2_HGCal_layer_"+std::to_string(layer)).c_str(), &impactsY_associatedChi2[layer]);
        } else if (_extrapolationType == DWC) {
            tree_->Branch(("impactX_HGCal_layer_"+std::to_string(layer)).c_str(), &impactX[layer]);
            tree_->Branch(("impactY_HGCal_layer_"+std::to_string(layer)).c_str(), &impactY[layer]);
        }
    }  


    if (_extrapolationType == DATURA) {
        tree_->Branch("kinkAngleX_DUT1", &kinkAngleX_DUT1);  
        tree_->Branch("kinkAngleY_DUT1", &kinkAngleY_DUT1);  
    } else if (_extrapolationType == DWC) {
        tree_->Branch("trackChi2_X", &impactX_associatedChi2[0]);  
        tree_->Branch("trackChi2_Y", &impactY_associatedChi2[0]);  
        tree_->Branch("dwcReferenceType", &dwcReferenceType);  
        tree_->Branch("m_x", &m_x);  
        tree_->Branch("m_y", &m_y);  
        tree_->Branch("b_x", &b_x);  
        tree_->Branch("b_y", &b_y);  
    }
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

    if (_extrapolationType == DATURA) {
        edm::Handle<std::vector<HGCalTBDATURATelescopeData> > daturatracks; 
        event.getByToken(DATURATrackToken, daturatracks);
        nTrackCounter=0;
        
        for(auto daturatrack : *daturatracks) {
            nTrackCounter++;
            for(int layer=1; layer<=m_nLayers; layer++) {        
                impactsX[layer].push_back(daturatrack.Extrapolation_XY(layer).first);
                impactsY[layer].push_back(daturatrack.Extrapolation_XY(layer).second);
                impactsX_associatedChi2[layer].push_back(daturatrack.Extrapolation_XY_Chi2(layer).first);
                impactsY_associatedChi2[layer].push_back(daturatrack.Extrapolation_XY_Chi2(layer).second);
            } 
            if(daturatrack.floatUserRecords.has("kinkAngleX_DUT1")) kinkAngleX_DUT1.push_back(daturatrack.floatUserRecords.get("kinkAngleX_DUT1"));
            if(daturatrack.floatUserRecords.has("kinkAngleY_DUT1")) kinkAngleY_DUT1.push_back(daturatrack.floatUserRecords.get("kinkAngleY_DUT1"));

        }
    } else if (_extrapolationType == DWC) {
        edm::Handle<HGCalTBDWCTrack> dwctrack;
        event.getByToken(DWCTrackToken, dwctrack);
        if (! dwctrack->valid) {
            nTrackCounter=0;
            dwcReferenceType=0;
            m_x = m_y = -999;
            b_x = b_y = -999;
            impactX_associatedChi2[0] = -999;
            impactY_associatedChi2[0] = -999;
            for(int layer=1; layer<=m_nLayers; layer++) {        
                impactX[layer] = -999;
                impactY[layer] = -999;
            }               
        }
        else {
            nTrackCounter=1;
            dwcReferenceType = dwctrack->referenceType;
            m_x = dwctrack->m_x;
            m_y = dwctrack->m_y;
            b_x = dwctrack->b_x;
            b_y = dwctrack->b_y;
            impactX_associatedChi2[0] = dwctrack->chi2_x;
            impactY_associatedChi2[0] = dwctrack->chi2_y;
            for(int layer=1; layer<=m_nLayers; layer++) {        
                impactX[layer] = (dwctrack->DWCExtrapolation_XY(layer).first);
                impactY[layer] = (dwctrack->DWCExtrapolation_XY(layer).second);
            }   

        } 
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
