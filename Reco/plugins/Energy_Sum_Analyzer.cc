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
#include "TH2Poly.h"
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

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "TFile.h"
#include "TTree.h"
  
//configuration1 based on Shilpis plot from 26.07.17: 
//https://indico.cern.ch/event/656159/contributions/2674177/attachments/1499439/2334623/update_simulation_geom.pdf
// 2 layers in EE (first one in the graphics was removed) and 4 in FH, indications are in mm

#define MAXVERTICES 6
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
		void InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV);


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
		std::vector<int> cellID_mostIntense_layer;
		std::vector<double> energyAll_layer, energyE1_layer, energyE7_layer, energyE19_layer;
		
		//quantitites for investigator wit event displays
		HGCalTBTopology IsCellValid;
		HGCalTBCellVertices TheCell;
		std::vector<std::pair<double, double>> CellXY;
		std::pair<double, double> CellCentreXY;

		std::vector<TH2Poly*> h_RecHit_MostIntense_layer;
		std::vector<TH2Poly*> h_RecHit_First7_layer;
		std::vector<TH2Poly*> h_RecHit_First19_layer;
		std::vector<TH2Poly*> h_RecHit_Occupancy_layer;
		std::vector<TH2Poly*> h_RecHit_Energy_layer;
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
		cellID_mostIntense_layer.push_back(0);
		energyAll_layer.push_back(0.);
		energyE1_layer.push_back(0.);
		energyE7_layer.push_back(0.);
		energyE19_layer.push_back(0.);
	}
	for (int l=1; l<=nLayers; l++) {
		outTree->Branch(("cellID_mostIntense_layer"+std::to_string(l)).c_str(), &cellID_mostIntense_layer[l-1], ("cellID_mostIntense_layer"+std::to_string(l)+"/I").c_str());
		outTree->Branch(("energyAll_layer"+std::to_string(l)).c_str(), &energyAll_layer[l-1], ("energyAll_layer"+std::to_string(l)+"/D").c_str());
		outTree->Branch(("energyE1_layer"+std::to_string(l)).c_str(), &energyE1_layer[l-1], ("energyE1_layer"+std::to_string(l)+"/D").c_str());
		outTree->Branch(("energyE7_layer"+std::to_string(l)).c_str(), &energyE7_layer[l-1], ("energyE7_layer"+std::to_string(l)+"/D").c_str());
		outTree->Branch(("energyE19_layer"+std::to_string(l)).c_str(), &energyE19_layer[l-1], ("energyE19_layer"+std::to_string(l)+"/D").c_str());
	}
	
	std::ostringstream os( std::ostringstream::ate );

	for( int ilayer=0; ilayer<nLayers; ilayer++ ){
		os.str("");
		os << "MostIntense_Layer" << (ilayer+1);
		h_RecHit_MostIntense_layer.push_back(fs->make<TH2Poly>());
		h_RecHit_MostIntense_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_MostIntense_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_MostIntense_layer[ilayer],ilayer,0,0);

		os.str("");
		os << "First7_Layer" << (ilayer+1);
		h_RecHit_First7_layer.push_back(fs->make<TH2Poly>());
		h_RecHit_First7_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_First7_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_First7_layer[ilayer],ilayer,0,0);

		os.str("");
		os << "First19_Layer" << (ilayer+1);
		h_RecHit_First19_layer.push_back(fs->make<TH2Poly>());
		h_RecHit_First19_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_First19_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_First19_layer[ilayer],ilayer,0,0);   

		os.str("");
		os << "Occupancy_Layer" << (ilayer+1);
		h_RecHit_Occupancy_layer.push_back(fs->make<TH2Poly>());
		h_RecHit_Occupancy_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_Occupancy_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_Occupancy_layer[ilayer],ilayer,0,0);   

		os.str("");
		os << "Energy_Layer" << (ilayer+1);
		h_RecHit_Energy_layer.push_back(fs->make<TH2Poly>());
		h_RecHit_Energy_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_Energy_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_Energy_layer[ilayer],ilayer,0,0);  
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
		int skiroc = Rechit.skiroc()+(layer-1)*4;  
		double iux = Rechit.getCellCenterCartesianCoordinate(0);
		double ivy = Rechit.getCellCenterCartesianCoordinate(1);
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap(layer);
			Sensors[layer]->setSensorSize(SensorSize);
		}

//		Sensors[layer]->addHit(Rechit, ADC_per_MIP[skiroc]);		//with MIP calibration

		if (Rechit.energyHigh() > MIP_cut)	{//only add if energy is higher than the MIP cut
			Sensors[layer]->addHit(Rechit, ADC_per_MIP[skiroc]);		//without MIP calibration
			h_RecHit_Occupancy_layer[layer-1]->Fill(iux, ivy);
			h_RecHit_Energy_layer[layer-1]->Fill(iux, ivy, Rechit.energyLow());
		}
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
	std::vector<std::pair<double, double> > relevantHitPositions;

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
		cellID_mostIntense_layer[it->first-1] = it->second->getMostIntensiveHit();
		energyE1_tot += it->second->getTotalWeight();
		energyE1_layer[it->first-1] = it->second->getTotalWeight();

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		for (size_t i=0; i<relevantHitPositions.size(); i++) {
			h_RecHit_MostIntense_layer[it->first-1]->Fill(relevantHitPositions[i].first/10., relevantHitPositions[i].second/10.);
		}
		relevantHitPositions.clear();
	}

	//step 4: sum of cells with highest intensity + one ring in a layer
	energyE7_tot = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(CONSIDERSEVEN, LINEARWEIGHTING);
		energyE7_tot += it->second->getTotalWeight();
		energyE7_layer[it->first-1] = it->second->getTotalWeight();

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		for (size_t i=0; i<relevantHitPositions.size(); i++) {
			h_RecHit_First7_layer[it->first-1]->Fill(relevantHitPositions[i].first/10., relevantHitPositions[i].second/10.);
		}
		relevantHitPositions.clear();
	}

	//step 5: sum of cells with highest intensity + two rings in a layer
	energyE19_tot = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(CONSIDERNINETEEN, LINEARWEIGHTING);
		energyE19_tot += it->second->getTotalWeight();
		energyE19_layer[it->first-1] = it->second->getTotalWeight();

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		for (size_t i=0; i<relevantHitPositions.size(); i++) {
			h_RecHit_First19_layer[it->first-1]->Fill(relevantHitPositions[i].first/10., relevantHitPositions[i].second/10.);
		}
		relevantHitPositions.clear();
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

void Energy_Sum_Analyzer::InitTH2Poly(TH2Poly& poly, int layerID, int sensorIU, int sensorIV)
{
  double HexX[MAXVERTICES] = {0.};
  double HexY[MAXVERTICES] = {0.};

  for(int iv = -7; iv < 8; iv++) {
    for(int iu = -7; iu < 8; iu++) {
      if(!IsCellValid.iu_iv_valid(layerID, sensorIU, sensorIV, iu, iv, 128)) continue;
      CellXY = TheCell.GetCellCoordinatesForPlots(layerID, sensorIU, sensorIV, iu, iv, 128);
      assert(CellXY.size() == 4 || CellXY.size() == 6);
      unsigned int iVertex = 0;
      for(std::vector<std::pair<double, double>>::const_iterator it = CellXY.begin(); it != CellXY.end(); it++) {
	HexX[iVertex] =  it->first;
	HexY[iVertex] =  it->second;
	++iVertex;
      }
      //Somehow cloning of the TH2Poly was not working. Need to look at it. Currently physically booked another one.
      poly.AddBin(CellXY.size(), HexX, HexY);
    }//loop over iu
  }//loop over iv
}

void Energy_Sum_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Energy_Sum_Analyzer);