/* 
 * A module that writes a ROOT tuple containing only relevant information for imaging
 * of rechHits.
 */

/**
	@Author: Thorben Quast <tquast>
		24 Febr 2017
		thorben.quast@cern.ch / thorben.quast@rwth-aachen.de
*/


// system include files
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"

#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TTree.h"

                       
class Imaging_Tuple_Writer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit Imaging_Tuple_Writer(const edm::ParameterSet&);
		~Imaging_Tuple_Writer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;


		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;
		edm::EDGetTokenT<RunData> RunDataToken;	
		
		struct {
  			HGCalElectronicsMap emap_;
		} essource_;
		std::string _e_mapFile;


		std::vector<double> ADC_per_MIP;

		int eventCounter, run;

		//stuff to be written to the tree
		TTree* outTree;
		double pBeam;
		std::vector<uint32_t> layers;
		std::vector<uint32_t> types;
		std::vector<int32_t> ius;
		std::vector<int32_t> ivs;
		std::vector<double> xs;
		std::vector<double> ys;
		std::vector<double> energies;
};

Imaging_Tuple_Writer::Imaging_Tuple_Writer(const edm::ParameterSet& iConfig) {	
	
	// initialization
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");
	
	_e_mapFile = iConfig.getParameter<std::string>("e_mapFile_CERN");	
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_e_mapFile);
  	if (!io.load(fip.fullPath(), essource_.emap_)) {
	  throw cms::Exception("Unable to load electronics map");
	};

	eventCounter = 0;

	//initialize tree and set Branch addresses
	outTree = fs->make<TTree>("recHits", "recHits");
	outTree->Branch("pBeam", &pBeam, "pBeam/D");
	outTree->Branch("layers", &layers);
	outTree->Branch("types", &types);
	outTree->Branch("ius", &ius);
	outTree->Branch("ivs", &ivs);
	outTree->Branch("xs", &xs);
	outTree->Branch("ys", &ys);
	outTree->Branch("energies", &energies);

}//constructor ends here

Imaging_Tuple_Writer::~Imaging_Tuple_Writer() {
	return;
}

// ------------ method called for each event  ------------
void Imaging_Tuple_Writer::analyze(const edm::Event& event, const edm::EventSetup& setup) {

	edm::Handle<RunData> rd;
 	//get the relevant event information
	event.getByToken(RunDataToken, rd);

	int evId = event.id().event();
	run = rd->run;
	pBeam = rd->trueEnergy;
	if (rd->hasDanger) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" is skipped because it has DANGER=true"<<std::endl;
		return;
	}
	

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);


	for(auto Rechit : *Rechits) {	
		layers.push_back((Rechit.id()).layer());
		types.push_back((Rechit.id()).cellType());
		ius.push_back((Rechit.id()).iu());
		ivs.push_back((Rechit.id()).iv());
		xs.push_back(Rechit.getCellCenterCartesianCoordinate(0));
		ys.push_back(Rechit.getCellCenterCartesianCoordinate(1));


		uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
		HGCalTBElectronicsId eid(EID);	 
		int skiRocIndex = (eid.iskiroc() - 1) > 0 ? eid.iskiroc() - 1 : 0;	
		energies.push_back(Rechit.energy()/ADC_per_MIP[skiRocIndex]);
	}

	outTree->Fill();

	layers.clear();
	types.clear();
	ius.clear();
	ivs.clear();
	xs.clear();
	ys.clear();
	energies.clear();


}// analyze ends here

void Imaging_Tuple_Writer::beginJob() {	
}

void Imaging_Tuple_Writer::endJob() {

}

void Imaging_Tuple_Writer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Imaging_Tuple_Writer);