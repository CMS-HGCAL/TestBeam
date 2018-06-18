#include <iostream>
#include "TTree.h"
#include <fstream>
#include <sstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h" //for the runData type definition

#include <iomanip>
#include <set>

class MetaDataNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit MetaDataNtupler(const edm::ParameterSet&);
    ~MetaDataNtupler();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void endJob() override;
    void analyze(const edm::Event& , const edm::EventSetup&) override;

    // conditions
    std::string m_metaDataFile;
    
    edm::EDGetTokenT<RunData> RunDataToken; 
    
    
    std::map<int, std::vector<float> > metaData;

    // Output tree
    TTree* tree_;

    void clearVariables(); // function to clear tree variables/vectors



    // Variables for branches

    // event info
    unsigned int ev_run_;
    unsigned int ev_event_;

    int pdgID;
    double beamEnergy;
    int configuration;

    float biasCurrentCh0;
    float biasCurrentCh1;
    float biasCurrentCh2;
    float biasCurrentCh3;
    float humidity_RHDP4;
    float TCassette07;
    float tablePositionY;
    float humidity_air;
    float temperature_air;



};

void MetaDataNtupler::clearVariables(){
    // event info
    ev_run_ = 0;
    ev_event_ = 0;
    pdgID = 0;
    beamEnergy = 0;
    configuration = 0;

    biasCurrentCh0 = 0;
    biasCurrentCh1 = 0;
    biasCurrentCh2 = 0;
    biasCurrentCh3 = 0;
    humidity_RHDP4 = 0;
    TCassette07 = 0;
    tablePositionY = 0;
    humidity_air = 0;
    temperature_air = 0;    
};

MetaDataNtupler::MetaDataNtupler(const edm::ParameterSet& iConfig) :
    m_metaDataFile(iConfig.getUntrackedParameter<std::string>("metaDataFile",""))
{
    RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
    

    usesResource("TFileService");
    edm::Service<TFileService> fs;    
    tree_ = fs->make<TTree>("meta", "meta");

    // event info
    tree_->Branch("event", &ev_event_);
    tree_->Branch("run", &ev_run_);

    tree_->Branch("pdgID", &pdgID);
    tree_->Branch("beamEnergy", &beamEnergy);

    tree_->Branch("configuration", &configuration);

    tree_->Branch("biasCurrentCh0", &biasCurrentCh0);
    tree_->Branch("biasCurrentCh1", &biasCurrentCh1);
    tree_->Branch("biasCurrentCh2", &biasCurrentCh2);
    tree_->Branch("biasCurrentCh3", &biasCurrentCh3);
    tree_->Branch("humidity_RHDP4", &humidity_RHDP4);
    tree_->Branch("TCassette07", &TCassette07);
    tree_->Branch("tablePositionY", &tablePositionY);
    tree_->Branch("humidity_air", &humidity_air);
    tree_->Branch("temperature_air", &temperature_air);

    //reads the meta data file
    std::fstream file; 
    char fragment[100];
    int readCounter = -1;

    readCounter = -1;
    file.open(m_metaDataFile.c_str(), std::fstream::in);
    std::cout<<"Reading file "<<m_metaDataFile<<" -open: "<<file.is_open()<<std::endl;
    int run=0;
    
    while (file.is_open() && !file.eof()) {
        readCounter++;
        file >> fragment;
        if (readCounter==0) run=atoi(fragment);
        else if (readCounter==9) {
            metaData[run].push_back(atof(fragment));
            readCounter=-1;
        } else {
            metaData[run].push_back(atof(fragment));
        }
    }

}


MetaDataNtupler::~MetaDataNtupler()
{

}

void MetaDataNtupler::beginJob()
{
}

void MetaDataNtupler::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    clearVariables();

    edm::Handle<RunData> rd;
    event.getByToken(RunDataToken, rd);   

    ev_run_ = rd->run;
    ev_event_ = rd->event;

    pdgID = rd->pdgID;
    beamEnergy = rd->energy;
    configuration = rd->configuration;

    biasCurrentCh0 = metaData[ev_run_][0];
    biasCurrentCh1 = metaData[ev_run_][1];
    biasCurrentCh2 = metaData[ev_run_][2];
    biasCurrentCh3 = metaData[ev_run_][3];
    humidity_RHDP4 = metaData[ev_run_][4];
    TCassette07 = metaData[ev_run_][5];
    tablePositionY = metaData[ev_run_][6];
    humidity_air = metaData[ev_run_][7];
    temperature_air = metaData[ev_run_][8];

    tree_->Fill();
}


void MetaDataNtupler::endJob()
{
}

void MetaDataNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MetaDataNtupler);
