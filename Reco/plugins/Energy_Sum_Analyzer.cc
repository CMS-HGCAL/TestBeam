/* 
 * Determination of the position resolution of the setup.
 */

/**
	@Author: Thorben Quast <tquast>
		18 September 2017
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
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "TFile.h"
#include "TTree.h"
  
//configuration1 based on Shilpis plot from 26.07.17: 
//https://indico.cern.ch/event/656159/contributions/2674177/attachments/1499439/2334623/update_simulation_geom.pdf
// 2 layers in EE (first one in the graphics was removed) and 4 in FH, indications are in mm


//#define DEBUG

class Energy_Sum_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit Energy_Sum_Analyzer(const edm::ParameterSet&);
		~Energy_Sum_Analyzer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;


		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;
	 	
		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<WireChambers> MWCToken;
		
				
		std::vector<double> ADC_per_MIP;		//one value per skiroc
		int LayersConfig;
		int SensorSize;
		int nLayers;

		double MIP_cut;

		//helper variables that are set within the event loop, i.e. are defined per event
		std::map<int, SensorHitMap*> Sensors;

		//stuff to be written to the tree
		TTree* outTree;
		int configuration, evId, eventCounter, run, pdgID; 	//eventCounter: counts the events in this analysis run to match information within ove event to each other
		double beamEnergy;		//energy of the beam particle

		double energyAll_tot, energyE1_tot, energyE7_tot, energyE19_tot;
		std::vector<double> energyAll_layer, energyE1_layer, energyE7_layer, energyE19_layer;
		

};

Energy_Sum_Analyzer::Energy_Sum_Analyzer(const edm::ParameterSet& iConfig) {	
	
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	
	//read the layer configuration
	LayersConfig = iConfig.getParameter<int>("layers_config");


	eventCounter = 0;

	SensorSize = iConfig.getParameter<int>("SensorSize");
	nLayers = iConfig.getParameter<int>("nLayers");
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");

	MIP_cut = 4.;		//todo: make as configuration
	MIP_cut *= 49.3;	//50 HG ADC counts are one MIP --> cell must at least have that value in order to be summed.

	//initialize tree and set Branch addresses
	outTree = fs->make<TTree>("energySums", "energySums");
	outTree->Branch("configuration", &configuration, "configuration/I");
	outTree->Branch("eventId", &evId, "eventId/I");	//event ID as it comes from the reader, as it is stored in the txt files
	outTree->Branch("eventCounter", &eventCounter, "eventCounter/I");	//event counter, current iteration, indexing occurs chronologically in the readin plugins
	outTree->Branch("run", &run, "run/I");
	outTree->Branch("pdgID", &pdgID, "pdgID/I");
	outTree->Branch("beamEnergy", &beamEnergy, "beamEnergy/D");	//electron energy in GeV
	outTree->Branch("energyAll_tot", &energyAll_tot, "energyAll_tot/D");
	outTree->Branch("energyE1_tot", &energyE1_tot, "energyE1_tot/D");
	outTree->Branch("energyE7_tot", &energyE7_tot, "energyE7_tot/D");
	outTree->Branch("energyE19_tot", &energyE19_tot, "energyE19_tot/D");

	for (int l=1; l<=nLayers; l++) {
		energyAll_layer.push_back(0.);
		energyE1_layer.push_back(0.);
		energyE7_layer.push_back(0.);
		energyE19_layer.push_back(0.);
	}
	for (int l=1; l<=nLayers; l++) {
		outTree->Branch(("energyAll_layer"+std::to_string(l)).c_str(), &energyAll_layer[l-1], ("energyAll_layer"+std::to_string(l)+"/D").c_str());
		outTree->Branch(("energyE1_layer"+std::to_string(l)).c_str(), &energyE1_layer[l-1], ("energyE1_layer"+std::to_string(l)+"/D").c_str());
		outTree->Branch(("energyE7_layer"+std::to_string(l)).c_str(), &energyE7_layer[l-1], ("energyE7_layer"+std::to_string(l)+"/D").c_str());
		outTree->Branch(("energyE19_layer"+std::to_string(l)).c_str(), &energyE19_layer[l-1], ("energyE19_layer"+std::to_string(l)+"/D").c_str());
	}
	

}//constructor ends here

Energy_Sum_Analyzer::~Energy_Sum_Analyzer() {
	return;
}

// ------------ method called for each event  ------------
void Energy_Sum_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {

	edm::Handle<RunData> rd;
 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	configuration = rd->configuration;
	evId = event.id().event();
	run = rd->run;
	pdgID = std::atoi( (rd->runType).c_str() );
	eventCounter = rd->event;
	beamEnergy = rd->energy;
	
	if (rd->booleanUserRecords.has("hasDanger")&&rd->booleanUserRecords.get("hasDanger")) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" ("<<beamEnergy<<"GeV)  is skipped because somthing went wrong"<<std::endl;
		return;
	}

	if (run == -1) {
		std::cout<<"Run is not in configuration file - is ignored."<<std::endl;
		return;
	}


	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);


	//step 1: Reduce the information to energy deposits/hits in x,y per sensor/layer 
	//fill the rechits:
	for(auto Rechit : *Rechits) {	
		int layer = (Rechit.id()).layer();
		//int skiroc = Rechit.skiroc()+(layer-1)*4;  
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap(layer);
			Sensors[layer]->setSensorSize(SensorSize);
		}

//		Sensors[layer]->addHit(Rechit, ADC_per_MIP[skiroc]);		//with MIP calibration
		
		if (Rechit.energy() > MIP_cut)	//only add if energy is higher than the MIP cut
			Sensors[layer]->addHit(Rechit, 1.);		//without MIP calibration
	}

	#ifdef DEBUG
		std::cout<<"run: "<<rd->run<<"  energy: "<<rd->energy<<"  type:" << rd->runType<<"   eventCounter: "<<rd->event<<std::endl;
	#endif
		
	
	//Possible event selection: sum of energies of all cells(=hits) from RecHits Collection and Clusters
	//sumEnergy = 0.;
	//for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {	
	//	sumEnergy += it->second->getTotalEnergy();
	//}

	//	considerationMethod = CONSIDERALL;
	//	considerationMethod = CONSIDERSEVEN;
	//	considerationMethod = CONSIDERNINETEEN;
	//	weightingMethod = LINEARWEIGHTING;
	//	weightingMethod = MOSTINTENSIVE;

	//step 2: sum of all cells in a layer
	energyAll_tot = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(CONSIDERALL, LINEARWEIGHTING);
		energyAll_tot += it->second->getTotalWeight();
		energyAll_layer[it->first-1] = it->second->getTotalWeight();
	}

	//step 3: sum of cells with highest intensity in a layer
	energyE1_tot = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(CONSIDERALL, MOSTINTENSIVE);
		energyE1_tot += it->second->getTotalWeight();
		energyE1_layer[it->first-1] = it->second->getTotalWeight();
	}

	//step 4: sum of cells with highest intensity + one ring in a layer
	energyE7_tot = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(CONSIDERSEVEN, LINEARWEIGHTING);
		energyE7_tot += it->second->getTotalWeight();
		energyE7_layer[it->first-1] = it->second->getTotalWeight();
	}

	//step 5: sum of cells with highest intensity + two rings in a layer
	energyE19_tot = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(CONSIDERNINETEEN, LINEARWEIGHTING);
		energyE19_tot += it->second->getTotalWeight();
		energyE19_layer[it->first-1] = it->second->getTotalWeight();
		#ifdef DEBUG
			std::cout<<"Roll position: "<<(it->first-1)<<" energyE19_layer = "<<energyE19_layer[it->first-1]<<std::endl;
		#endif
	}

	outTree->Fill();
	//fill the tree


	
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();



	
}// analyze ends here

void Energy_Sum_Analyzer::beginJob() {	
}

void Energy_Sum_Analyzer::endJob() {
	
}

void Energy_Sum_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Energy_Sum_Analyzer);