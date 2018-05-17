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

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include <iomanip>
#include <set>

class RecHitNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit RecHitNtupler(const edm::ParameterSet&);
    ~RecHitNtupler();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void endJob() override;
    void analyze(const edm::Event& , const edm::EventSetup&) override;

    // conditions
    std::string m_electronicMap;
    std::string m_detectorLayoutFile;

    struct {
	HGCalElectronicsMap emap_;
	HGCalTBDetectorLayout layout_;
    } essource_;

    // parameters
    int m_sensorsize;
    bool m_eventPlotter;
    int m_evtID;
    double m_mipThreshold;
    double m_noiseThreshold;

    // ---------- member data ---------------------------
    edm::EDGetTokenT<HGCalTBRecHitCollection> m_HGCalTBRecHitCollection;

    HGCalTBTopology IsCellValid;
    HGCalTBCellVertices TheCell;
    std::vector<std::pair<double, double>> CellXY;
    std::pair<double, double> CellCentreXY;

    // Output tree
    TTree* tree_;

    void clearVariables(); // function to clear tree variables/vectors

    // Variables for branches

    // event info
    unsigned int ev_run_;
    unsigned int ev_event_;

    // rechits
    std::vector<unsigned int> rechit_detid_;
    std::vector<unsigned int> rechit_module_;
    std::vector<unsigned int> rechit_layer_;
    std::vector<float> rechit_x_;
    std::vector<float> rechit_y_;
    std::vector<float> rechit_z_;
    std::vector<int> rechit_iu_;
    std::vector<int> rechit_iv_;
    std::vector<float> rechit_energy_;
    std::vector<float> rechit_energyHigh_;
    std::vector<float> rechit_energyLow_;
    std::vector<float> rechit_energyTot_;
    std::vector<float> rechit_time_;

};

void RecHitNtupler::clearVariables(){
    // event info
    ev_run_ = 0;
    ev_event_ = 0;

    // rechits
    rechit_detid_.clear();
    rechit_module_.clear();
    rechit_layer_.clear();
    rechit_x_.clear();
    rechit_y_.clear();
    rechit_z_.clear();
    rechit_iu_.clear();
    rechit_iv_.clear();
    rechit_energy_.clear();
    rechit_energyHigh_.clear();
    rechit_energyLow_.clear();
    rechit_energyTot_.clear();
    rechit_time_.clear();

};

RecHitNtupler::RecHitNtupler(const edm::ParameterSet& iConfig) :
    m_electronicMap(iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt")),
    m_detectorLayoutFile(iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt")),
    m_sensorsize(iConfig.getUntrackedParameter<int>("SensorSize",128)),
    m_eventPlotter(iConfig.getUntrackedParameter<bool>("EventPlotter",false)),
    m_mipThreshold(iConfig.getUntrackedParameter<double>("MipThreshold",5.0)),
    m_noiseThreshold(iConfig.getUntrackedParameter<double>("NoiseThreshold",0.5))
{
    m_HGCalTBRecHitCollection = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("InputCollection"));

    m_evtID=0;

    std::cout << iConfig.dump() << std::endl;

    HGCalCondObjectTextIO io(0);
    edm::FileInPath fip(m_electronicMap);
    if (!io.load(fip.fullPath(), essource_.emap_)) {
	throw cms::Exception("Unable to load electronics map");
    };
    fip=edm::FileInPath(m_detectorLayoutFile);
    if (!io.load(fip.fullPath(), essource_.layout_)) {
	throw cms::Exception("Unable to load detector layout file");
    };

    usesResource("TFileService");
    edm::Service<TFileService> fs;

    // Define tree and branches
    tree_ = fs->make<TTree>("hits", "HGC rechits");

    // event info
    tree_->Branch("event", &ev_event_);
    tree_->Branch("event", &ev_run_);

    // rechit
    tree_->Branch("rechit_detid",&rechit_detid_);
    tree_->Branch("rechit_module",&rechit_module_);
    tree_->Branch("rechit_layer",&rechit_layer_);
    tree_->Branch("rechit_x",&rechit_x_);
    tree_->Branch("rechit_y",&rechit_y_);
    tree_->Branch("rechit_z",&rechit_z_);
    tree_->Branch("rechit_iu",&rechit_iu_);
    tree_->Branch("rechit_iv",&rechit_iv_);
    tree_->Branch("rechit_energy",&rechit_energy_);
    tree_->Branch("rechit_energyHigh",&rechit_energyHigh_);
    tree_->Branch("rechit_energyLow",&rechit_energyLow_);
    tree_->Branch("rechit_energyTot",&rechit_energyTot_);
    tree_->Branch("rechit_time",&rechit_time_);

}


RecHitNtupler::~RecHitNtupler()
{

}

void RecHitNtupler::beginJob()
{
}

void RecHitNtupler::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
    clearVariables();

    usesResource("TFileService");
    edm::Service<TFileService> fs;

    edm::Handle<HGCalTBRecHitCollection> rhits;
    event.getByToken(m_HGCalTBRecHitCollection, rhits);

    ev_run_ = event.id().run();
    ev_event_ = event.id().event();

    for( auto hit : *rhits ){

	// get electronics channel
	HGCalTBElectronicsId eid( essource_.emap_.detId2eid( hit.id().rawId() ) );

	// get geometric channel
	if ( !IsCellValid.iu_iv_valid(
		 hit.id().layer(),
		 hit.id().sensorIU(),
		 hit.id().sensorIV(),
		 hit.id().iu(),hit.id().iv(),m_sensorsize)
	    ) continue;

	CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(
	    hit.id().layer(),
	    hit.id().sensorIU(),
	    hit.id().sensorIV(),
	    hit.id().iu(),hit.id().iv(),m_sensorsize
	    );

	// get layer and module ID
	HGCalTBLayer layer = essource_.layout_.at(hit.id().layer()-1);
	int moduleId = layer.at( hit.id().sensorIU(),hit.id().sensorIV() ).moduleID();

	// Fill hit info and position
	rechit_detid_.push_back(hit.id());
	// rechit_chip_ = eid.iskiroc();
	// rechit_channel_ = eid.ichannel();
	rechit_module_.push_back(moduleId);
	rechit_layer_.push_back(hit.id().layer());

	rechit_x_.push_back( CellCentreXY.first );
	rechit_y_.push_back( CellCentreXY.second );

	/*
	// or instead?
	rechit_x_.push_back( hit.cellCenter_x );
	rechit_y_.push_back( hit.cellCenter_y );
	*/

	// not available for now
	//rechit_z_.push_back();

	rechit_iu_.push_back( hit.id().iu() );
	rechit_iv_.push_back( hit.id().iv() );

	// Hit energy and time
	rechit_energy_.push_back( hit.energy() );

	rechit_energyHigh_.push_back( hit.energyHigh() );
	rechit_energyLow_.push_back( hit.energyLow() );
	rechit_energyTot_.push_back( hit.energyTot() );

	rechit_time_.push_back( hit.time() );

    } // end rechit loop

    tree_->Fill();
}


void RecHitNtupler::endJob()
{
}

void RecHitNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecHitNtupler);
