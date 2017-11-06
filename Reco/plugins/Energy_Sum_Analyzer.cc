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
		int SensorSize;
		int nLayers;

		double MIP_cut;

		//helper variables that are set within the event loop, i.e. are defined per event
		std::map<int, SensorHitMap*> Sensors;

		//stuff to be written to the tree

		int configuration, evId, eventCounter, run, pdgID; 	//eventCounter: counts the events in this analysis run to match information within ove event to each other
		double beamEnergy;		//energy of the beam particle

		
		//quantitites for investigator wit event displays
		HGCalTBTopology IsCellValid;
		HGCalTBCellVertices TheCell;
		std::vector<std::pair<double, double>> CellXY;
		std::pair<double, double> CellCentreXY;


		double energyAll_tot, energyE1_tot, energyE7_tot, energyE19_tot;
		std::vector<int> cellID_mostIntense_layer;
		std::vector<double> energyAll_layer, energyE1_layer, energyE7_layer, energyE19_layer;
		std::vector<int> NAll_layer, NE1_layer, NE7_layer, NE19_layer;

		TH1F* h_energyAll_tot;
		TH1F* h_energyE1_tot;
		TH1F* h_energyE7_tot;
		TH1F* h_energyE19_tot;

		std::vector<TH1F*> h_energyAll_layer;
		std::vector<TH1F*> h_energyE1_layer;
		std::vector<TH1F*> h_energyE7_layer;
		std::vector<TH1F*> h_energyE19_layer;
		
		std::vector<TH1F*> h_NAll_layer;
		std::vector<TH1F*> h_NE1_layer;
		std::vector<TH1F*> h_NE7_layer;
		std::vector<TH1F*> h_NE19_layer;	

		std::vector<TH1F*> h_energyE1perE7_layer;
		std::vector<TH1F*> h_energyE1perE19_layer;
		std::vector<TH1F*> h_energyE1perAll_layer;
		std::vector<TH1F*> h_energyE7perE19_layer;
		std::vector<TH1F*> h_energyE7perAll_layer;
		std::vector<TH1F*> h_energyE19perAll_layer;
		
		std::vector<TH2Poly*> h_RecHit_MostIntense_layer;
		std::vector<TH2Poly*> h_RecHit_First7_layer;
		std::vector<TH2Poly*> h_RecHit_First19_layer;
		std::vector<TH2Poly*> h_RecHit_Occupancy_layer;
		std::vector<TH2Poly*> h_RecHit_Energy_layer;

		//for debugging data vs. simulation
		int cellIDToFocusOn;
};

Energy_Sum_Analyzer::Energy_Sum_Analyzer(const edm::ParameterSet& iConfig) {	
	cellIDToFocusOn = -1;

	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));

	eventCounter = 0;

	SensorSize = iConfig.getParameter<int>("SensorSize");
	nLayers = iConfig.getParameter<int>("nLayers");
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");

	MIP_cut = 4.;		
	
	std::ostringstream os( std::ostringstream::ate );

	os.str("");
	os << "EnergyAll";
	h_energyAll_tot = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 2500.*nLayers/6);
	h_energyAll_tot->SetName(os.str().c_str());
	h_energyAll_tot->SetTitle(os.str().c_str());
	os.str("");
	os << "EnergyE1";
	h_energyE1_tot = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 1000.*nLayers/6);
	h_energyE1_tot->SetName(os.str().c_str());
	h_energyE1_tot->SetTitle(os.str().c_str());
	os.str("");
	os << "EnergyE7";
	h_energyE7_tot = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 2000.*nLayers/6);
	h_energyE7_tot->SetName(os.str().c_str());
	h_energyE7_tot->SetTitle(os.str().c_str());
	os.str("");
	os << "EnergyE19";
	h_energyE19_tot = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 2500.*nLayers/6);
	h_energyE19_tot->SetName(os.str().c_str());
	h_energyE19_tot->SetTitle(os.str().c_str());

	
	for( int ilayer=0; ilayer<nLayers; ilayer++ ){
		cellID_mostIntense_layer.push_back(0);
		energyAll_layer.push_back(0.);
		energyE1_layer.push_back(0.);
		energyE7_layer.push_back(0.);
		energyE19_layer.push_back(0.);
		NAll_layer.push_back(0);
		NE1_layer.push_back(0);
		NE7_layer.push_back(0);
		NE19_layer.push_back(0);
		
		os.str("");
		os << "Layer" << (ilayer+1);
		TFileDirectory layer_dir = fs->mkdir( os.str().c_str() );

		os.str("");
		os << "EnergyE1_Layer" << (ilayer+1);
		h_energyE1_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 200.));
		h_energyE1_layer[ilayer]->SetName(os.str().c_str());
		h_energyE1_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "EnergyE7_Layer" << (ilayer+1);
		h_energyE7_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 400.));
		h_energyE7_layer[ilayer]->SetName(os.str().c_str());
		h_energyE7_layer[ilayer]->SetTitle(os.str().c_str());
		

		os.str("");
		os << "EnergyE19_Layer" << (ilayer+1);
		h_energyE19_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 500.));
		h_energyE19_layer[ilayer]->SetName(os.str().c_str());
		h_energyE19_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "EnergyAll_Layer" << (ilayer+1);
		h_energyAll_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 500.));
		h_energyAll_layer[ilayer]->SetName(os.str().c_str());
		h_energyAll_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "NE1_Layer" << (ilayer+1);
		h_NE1_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 2, -0.5, 1.5));
		h_NE1_layer[ilayer]->SetName(os.str().c_str());
		h_NE1_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "NE7_Layer" << (ilayer+1);
		h_NE7_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 8, -0.5, 7.5));
		h_NE7_layer[ilayer]->SetName(os.str().c_str());
		h_NE7_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "NE19_Layer" << (ilayer+1);
		h_NE19_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 20, -0.5, 19.5));
		h_NE19_layer[ilayer]->SetName(os.str().c_str());
		h_NE19_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "NAll_Layer" << (ilayer+1);
		h_NAll_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 31, -0.5, 30.5));
		h_NAll_layer[ilayer]->SetName(os.str().c_str());
		h_NAll_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "EnergyE1perE7_Layer" << (ilayer+1);
		h_energyE1perE7_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 1.));
		h_energyE1perE7_layer[ilayer]->SetName(os.str().c_str());
		h_energyE1perE7_layer[ilayer]->SetTitle(os.str().c_str());
		
		os.str("");
		os << "EnergyE1perE19_Layer" << (ilayer+1);
		h_energyE1perE19_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 1.));
		h_energyE1perE19_layer[ilayer]->SetName(os.str().c_str());
		h_energyE1perE19_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "EnergyE1perAll_Layer" << (ilayer+1);
		h_energyE1perAll_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 1.));
		h_energyE1perAll_layer[ilayer]->SetName(os.str().c_str());
		h_energyE1perAll_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "EnergyE7perE19_Layer" << (ilayer+1);
		h_energyE7perE19_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 1.));
		h_energyE7perE19_layer[ilayer]->SetName(os.str().c_str());
		h_energyE7perE19_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "EnergyE7perAll_Layer" << (ilayer+1);
		h_energyE7perAll_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 1.));
		h_energyE7perAll_layer[ilayer]->SetName(os.str().c_str());
		h_energyE7perAll_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "EnergyE19perAll_Layer" << (ilayer+1);
		h_energyE19perAll_layer.push_back(layer_dir.make<TH1F>(os.str().c_str(), os.str().c_str(), 50, 0., 1.));
		h_energyE19perAll_layer[ilayer]->SetName(os.str().c_str());
		h_energyE19perAll_layer[ilayer]->SetTitle(os.str().c_str());

		os.str("");
		os << "MostIntense_Layer" << (ilayer+1);
		h_RecHit_MostIntense_layer.push_back(layer_dir.make<TH2Poly>());
		h_RecHit_MostIntense_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_MostIntense_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_MostIntense_layer[ilayer],ilayer,0,0);

		os.str("");
		os << "First7_Layer" << (ilayer+1);
		h_RecHit_First7_layer.push_back(layer_dir.make<TH2Poly>());
		h_RecHit_First7_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_First7_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_First7_layer[ilayer],ilayer,0,0);

		os.str("");
		os << "First19_Layer" << (ilayer+1);
		h_RecHit_First19_layer.push_back(layer_dir.make<TH2Poly>());
		h_RecHit_First19_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_First19_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_First19_layer[ilayer],ilayer,0,0);   

		os.str("");
		os << "Occupancy_Layer" << (ilayer+1);
		h_RecHit_Occupancy_layer.push_back(layer_dir.make<TH2Poly>());
		h_RecHit_Occupancy_layer[ilayer]->SetName( os.str().c_str() );
		h_RecHit_Occupancy_layer[ilayer]->SetTitle( os.str().c_str() );
		InitTH2Poly(*h_RecHit_Occupancy_layer[ilayer],ilayer,0,0);   

		os.str("");
		os << "Energy_Layer" << (ilayer+1);
		h_RecHit_Energy_layer.push_back(layer_dir.make<TH2Poly>());
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
	pdgID = rd->pdgID;
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
		int skiroc = Rechit.skiroc();
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap(layer);
			Sensors[layer]->setSensorSize(SensorSize);
		}

		if (Rechit.energyHigh() > MIP_cut*ADC_per_MIP[skiroc])	{//only add if energy is higher than the MIP cut
			Sensors[layer]->addHit(Rechit, ADC_per_MIP[skiroc]);		//without MIP calibration
		}
	}

	#ifdef DEBUG
		std::cout<<"run: "<<rd->run<<"  energy: "<<rd->energy<<"  pdgID:" << pdgID<<"   eventCounter: "<<rd->event<<std::endl;
	#endif

	
	std::vector<std::pair<double, double> > relevantHitPositions;


	//step 2: sum of cells with highest intensity in a layer
	energyE1_tot = 0.;
	energyE7_tot = 0.;
	energyE19_tot = 0.;
	energyAll_tot = 0.;
	
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer

		it->second->calculateCenterPosition(CONSIDERALL, MOSTINTENSIVE);
		cellID_mostIntense_layer[it->first-1] = it->second->getMostIntensiveHit();
		energyE1_tot += it->second->getTotalWeight();
		energyE1_layer[it->first-1] = it->second->getTotalWeight();

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE1_layer[it->first-1] = (int)relevantHitPositions.size();
		if (cellIDToFocusOn==-1 || (cellID_mostIntense_layer[it->first-1] % 1000) == cellIDToFocusOn) {	
			for (size_t i=0; i<relevantHitPositions.size(); i++) {
				h_RecHit_MostIntense_layer[it->first-1]->Fill(relevantHitPositions[i].first/10., relevantHitPositions[i].second/10.);
			}
		}
		relevantHitPositions.clear();
	

		//step 3: sum of cells with highest intensity + one ring in a layer
		it->second->calculateCenterPosition(CONSIDERSEVEN, LINEARWEIGHTING);
		energyE7_tot += it->second->getTotalWeight();
		energyE7_layer[it->first-1] = it->second->getTotalWeight();

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE7_layer[it->first-1] = (int)relevantHitPositions.size();
		if (cellIDToFocusOn==-1 || (cellID_mostIntense_layer[it->first-1] % 1000) == cellIDToFocusOn) {
			for (size_t i=0; i<relevantHitPositions.size(); i++) {
				h_RecHit_First7_layer[it->first-1]->Fill(relevantHitPositions[i].first/10., relevantHitPositions[i].second/10.);
			}
		}
		relevantHitPositions.clear();
	

		//step 4: sum of cells with highest intensity + two rings in a layer
		it->second->calculateCenterPosition(CONSIDERNINETEEN, LINEARWEIGHTING);
		energyE19_tot += it->second->getTotalWeight();
		energyE19_layer[it->first-1] = it->second->getTotalWeight();

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NE19_layer[it->first-1] = (int)relevantHitPositions.size();
		if (cellIDToFocusOn==-1 || (cellID_mostIntense_layer[it->first-1] % 1000) == cellIDToFocusOn) {
			for (size_t i=0; i<relevantHitPositions.size(); i++) {
				h_RecHit_First19_layer[it->first-1]->Fill(relevantHitPositions[i].first/10., relevantHitPositions[i].second/10.);
			}
		}
		relevantHitPositions.clear();
	

		//step 5: sum of all cells in a layer
		it->second->calculateCenterPosition(CONSIDERALL, LINEARWEIGHTING);
		energyAll_tot += it->second->getTotalWeight();
		energyAll_layer[it->first-1] = it->second->getTotalWeight();

		relevantHitPositions = it->second->getHitPositionsForPositioning();
		NAll_layer[it->first-1] = (int)relevantHitPositions.size();
	}

	//fill the histograms here

	h_energyAll_tot->Fill(energyAll_tot);
	h_energyE1_tot->Fill(energyE1_tot);
	h_energyE7_tot->Fill(energyE7_tot);
	h_energyE19_tot->Fill(energyE19_tot);

	for (int l=0; l<nLayers; l++) {
		if (cellIDToFocusOn==-1 || (cellID_mostIntense_layer[l] % 1000) != 236) continue;
		h_energyAll_layer[l]->Fill(energyAll_layer[l]);
		h_energyE1_layer[l]->Fill(energyE1_layer[l]);
		h_energyE7_layer[l]->Fill(energyE7_layer[l]);
		h_energyE19_layer[l]->Fill(energyE19_layer[l]);
	
		h_NAll_layer[l]->Fill(NAll_layer[l]);
		h_NE1_layer[l]->Fill(NE1_layer[l]);
		h_NE7_layer[l]->Fill(NE7_layer[l]);
		h_NE19_layer[l]->Fill(NE19_layer[l]);


		if (NE1_layer[l]==0) continue;	//only fill ratios if here are rings to be computed
		h_energyE1perE7_layer[l]->Fill( energyE1_layer[l] / energyE7_layer[l]);
		h_energyE1perE19_layer[l]->Fill( energyE1_layer[l] / energyE19_layer[l]);
		h_energyE1perAll_layer[l]->Fill( energyE1_layer[l] / energyAll_layer[l]);
		h_energyE7perE19_layer[l]->Fill( energyE7_layer[l] / energyE19_layer[l]);
		h_energyE7perAll_layer[l]->Fill( energyE7_layer[l] / energyAll_layer[l]);
		h_energyE19perAll_layer[l]->Fill( energyE19_layer[l] / energyAll_layer[l]);
		
	}
	
	//fill occupancy plots
	for(auto Rechit : *Rechits) {	
		int layer = (Rechit.id()).layer();
		double iux = Rechit.getCellCenterCartesianCoordinate(0);
		double ivy = Rechit.getCellCenterCartesianCoordinate(1);
		int skiroc = Rechit.skiroc()+(layer-1)*4;  

		if ((cellID_mostIntense_layer[layer-1] % 1000) != 236) continue;
		if (Rechit.energyHigh() <= MIP_cut*ADC_per_MIP[skiroc]) continue;
		if ((Rechit.id()).cellType() != 0) continue;
		h_RecHit_Occupancy_layer[layer-1]->Fill(iux, ivy);
		double _energy = Rechit.energy() / ADC_per_MIP[skiroc];
  		if (Rechit.checkFlag(HGCalTBRecHit::kLowGainSaturated)) {
    		_energy = Rechit.energyLow() * 8. / ADC_per_MIP[skiroc] ;  
  		}
		h_RecHit_Energy_layer[layer-1]->Fill(iux, ivy, _energy);
	}


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