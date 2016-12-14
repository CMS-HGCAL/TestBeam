/* 
 * Determination of the position resolution of the setup.
 * Hits per layer are calculated using dedicated weighting algorithms.
 * Hypothetical particle tracks are obtained from all but one layer.
 * The predicted position in that layer is then compared to the reconstructed one.
 * Thus, for each event and each layer there is one deviation to be calculated and filled into a 2D histogram.
 */

/**
	@Author: Thorben Quast <tquast>
		16 Nov 2016
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
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TTree.h"
  
//configuration1:
double config1Positions[] = {0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92}; 	 //z-coordinate in cm
double config1X0Depths[] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42}; //in radiation lengths, copied from layerSumAnalyzer

//configuration2:
double config2Positions[] = {0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4}; 				//z-coordinate in cm
double config2X0Depths[] = {5.048, 3.412, 3.412, 2.866, 2.512, 1.625, 2.368, 6.021}; //in radiation lengths, copied from layerSumAnalyzer
                     
class Position_Resolution_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit Position_Resolution_Analyzer(const edm::ParameterSet&);
		~Position_Resolution_Analyzer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;
	 	edm::EDGetToken HGCalTBClusterCollection_Token;
  	edm::EDGetToken HGCalTBClusterCollection7_Token;
  	edm::EDGetToken HGCalTBClusterCollection19_Token;

		edm::EDGetTokenT<RunData> RunDataToken;	
		
		ConsiderationMethod considerationMethod;
		WeightingMethod weightingMethod;
		TrackFittingMethod fittingMethod;		
		FitPointWeightingMethod fitPointWeightingMethod;

		std::vector<int> EventsFor2DGraphs;
		double pedestalThreshold;
		std::vector<double> Layer_Z_Positions;
		std::vector<double> Layer_Z_X0s;
		std::vector<double> ADC_per_MIP;
		int LayersConfig;
		int nLayers;
		int SensorSize;

		double totalEnergyThreshold;

		int ClusterVetoCounter;
		int HitsVetoCounter;
		int CommonVetoCounter;

		std::map<int, int> successfulFitCounter, failedFitCounter;


		//helper variables that are set within the event loop, i.e. are defined per event
		std::map<int, SensorHitMap*> Sensors;
		std::map<int, ParticleTrack*> Tracks;
		std::vector<double> x_predicted_v,y_predicted_v, x_true_v, y_true_v, layerZ_cm_v;

		//stuff to be written to the tree
		TTree* outTree;
		int configuration, evId, run, layer;
		double energy;
		double layerWeight, sumFitWeights, sumEnergy, sumClusterEnergy;
		double x_predicted, x_predicted_err, y_predicted, y_predicted_err, x_true, x_true_err, y_true, y_true_err, deltaX, deltaY;
		double x_predicted_to_closest_cell, y_predicted_to_closest_cell, x_true_to_closest_cell, y_true_to_closest_cell, layerZ_cm, layerZ_X0, deviation;
};

Position_Resolution_Analyzer::Position_Resolution_Analyzer(const edm::ParameterSet& iConfig) {
	gStyle->SetOptStat();
	
	// initialization
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
  HGCalTBClusterCollection_Token = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS"));
  HGCalTBClusterCollection7_Token = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS7"));
  HGCalTBClusterCollection19_Token = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS19"));

	//read the cell consideration option to calculate the central hit point
	std::string methodString = iConfig.getParameter<std::string>("considerationMethod");
	if (methodString == "all")
  	considerationMethod = CONSIDERALL;
  else if (methodString == "closest7")
  	considerationMethod = CONSIDERSEVEN;
  else if (methodString == "closest19")
  	considerationMethod = CONSIDERNINETEEN;
  else if(methodString == "clustersAll")
  	considerationMethod = CONSIDERCLUSTERSALL;
  else if(methodString == "clusters7")
  	considerationMethod = CONSIDERCLUSTERSSEVEN;
  else if(methodString == "clusters19")
  	considerationMethod = CONSIDERCLUSTERSNINETEEN;

	//read the weighting method to obtain the central hit point
	methodString = iConfig.getParameter<std::string>("weightingMethod");
	if (methodString == "squaredWeighting")
		weightingMethod = SQUAREDWEIGHTING;	
	else if (methodString == "linearWeighting")
		weightingMethod = LINEARWEIGHTING;
	else if (methodString == "logWeighting_3.0_1.0")
		weightingMethod = LOGWEIGHTING_30_10;
	else if (methodString == "logWeighting_3.0_1.5")
		weightingMethod = LOGWEIGHTING_30_15;
	else if (methodString == "logWeighting_4.0_1.0")
		weightingMethod = LOGWEIGHTING_40_10;
	else if (methodString == "logWeighting_4.0_1.5")
		weightingMethod = LOGWEIGHTING_40_15;
	else if (methodString == "logWeighting_5.0_1.0")
		weightingMethod = LOGWEIGHTING_50_10;
	else if (methodString == "logWeighting_5.0_1.5")
		weightingMethod = LOGWEIGHTING_50_15;
	else if (methodString == "logWeighting_6.0_1.0")
		weightingMethod = LOGWEIGHTING_60_10;
	else if (methodString == "logWeighting_6.0_1.5")
		weightingMethod = LOGWEIGHTING_60_15;
	else if (methodString == "logWeighting_7.0_1.0")
		weightingMethod = LOGWEIGHTING_70_10;
	else if (methodString == "logWeighting_7.0_1.5")
		weightingMethod = LOGWEIGHTING_70_15;	
	else 
		weightingMethod = DEFAULTWEIGHTING;

	//read the track fitting method
	methodString = iConfig.getParameter<std::string>("fittingMethod");
	if (methodString == "lineTGraphErrors")
		fittingMethod = LINEFITTGRAPHERRORS;
	else if (methodString == "pol2TGraphErrors")
		fittingMethod = POL2TGRAPHERRORS;
	else if (methodString == "pol3TGraphErrors")
		fittingMethod = POL3TGRAPHERRORS;
	else 
		fittingMethod = DEFAULTFITTING;

	//read the fit point weighting technique:
	methodString = iConfig.getParameter<std::string>("fitPointWeightingMethod");
	if (methodString == "none")
		fitPointWeightingMethod = NONE;
	else if (methodString == "linear")
		fitPointWeightingMethod = LINEAR;
	else if (methodString == "squared")
		fitPointWeightingMethod = SQUARED;
	else if (methodString == "logarithmic")
		fitPointWeightingMethod = LOGARITHMIC;
	else if (methodString == "exponential")
		fitPointWeightingMethod = EXPONENTIAL;
	else 
		fitPointWeightingMethod = NONE;

	//read the layer configuration
	LayersConfig = iConfig.getParameter<int>("layers_config");
	if (LayersConfig == 1) {
		Layer_Z_Positions = std::vector<double>(config1Positions, config1Positions + sizeof(config1Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(config1X0Depths, config1X0Depths + sizeof(config1X0Depths)/sizeof(double));
	} if (LayersConfig == 2) {
		Layer_Z_Positions = std::vector<double>(config2Positions, config2Positions + sizeof(config2Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(config2X0Depths, config2X0Depths + sizeof(config2X0Depths)/sizeof(double));
	} else {
		Layer_Z_Positions = std::vector<double>(config1Positions, config1Positions + sizeof(config1Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(config1X0Depths, config1X0Depths + sizeof(config1X0Depths)/sizeof(double));
	}

	pedestalThreshold = iConfig.getParameter<double>("pedestalThreshold");
	SensorSize = iConfig.getParameter<int>("SensorSize");
	nLayers = iConfig.getParameter<int>("nLayers");
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");

	totalEnergyThreshold = iConfig.getParameter<double>("totalEnergyThreshold");

	//making 2DGraphs per event?
	EventsFor2DGraphs = iConfig.getParameter<std::vector<int> >("EventsFor2DGraphs");

	//initialize tree and set Branch addresses
	outTree = fs->make<TTree>("deviations", "deviations");
	outTree->Branch("configuration", &configuration, "configuration/I");
	outTree->Branch("eventId", &evId, "eventId/I");
	outTree->Branch("run", &run, "run/I");
	outTree->Branch("layer", &layer, "layer/I");
	outTree->Branch("energy", &energy, "energy/D");
	outTree->Branch("layerWeight", &layerWeight, "layerWeight/D");
	outTree->Branch("sumEnergy", &sumEnergy, "sumEnergy/D");
	outTree->Branch("sumClusterEnergy", &sumClusterEnergy, "sumClusterEnergy/D");
	outTree->Branch("sumFitWeights", &sumFitWeights, "sumFitWeights/D");
	outTree->Branch("x_predicted", &x_predicted, "x_predicted/D");
	outTree->Branch("x_predicted_to_closest_cell", &x_predicted_to_closest_cell, "x_predicted_to_closest_cell/D");
	outTree->Branch("x_predicted_err", &x_predicted_err, "x_predicted_err/D");
	outTree->Branch("y_predicted", &y_predicted, "y_predicted/D");
	outTree->Branch("y_predicted_to_closest_cell", &y_predicted_to_closest_cell, "y_predicted_to_closest_cell/D");
	outTree->Branch("y_predicted_err", &y_predicted_err, "y_predicted_err/D");
	outTree->Branch("x_true", &x_true, "x_true/D");
	outTree->Branch("x_true_to_closest_cell", &x_true_to_closest_cell, "x_true_to_closest_cell/D");
	outTree->Branch("x_true_err", &x_true_err, "x_true_err/D");
	outTree->Branch("y_true", &y_true, "y_true/D");
	outTree->Branch("y_true_to_closest_cell", &y_true_to_closest_cell, "y_true_to_closest_cell/D");
	outTree->Branch("y_true_err", &y_true_err, "y_true_err/D");
	outTree->Branch("deltaX", &deltaX, "deltaX/D");
	outTree->Branch("deltaY", &deltaY, "deltaY/D");
	outTree->Branch("layerZ_cm", &layerZ_cm, "layerZ_cm/D");
	outTree->Branch("layerZ_X0", &layerZ_X0, "layerZ_X0/D");
	outTree->Branch("deviation", &deviation, "deviation/D");


	ClusterVetoCounter = 0;
	HitsVetoCounter = 0;
	CommonVetoCounter = 0;
}//constructor ends here

Position_Resolution_Analyzer::~Position_Resolution_Analyzer() {
	return;
}

// ------------ method called for each event  ------------
void Position_Resolution_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
	edm::Handle<RunData> rd;

 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	configuration = rd->configuration;
	evId = event.id().event();
	run = rd->run;
	energy = rd->energy;
	if (run == -1) {
		std::cout<<"Run is not in configuration file - is ignored."<<std::endl;
		return;
	}

	//check if 2DGraphs are to be made
	bool make2DGraphs = false;
	std::vector<int>::iterator findPosition = std::find(EventsFor2DGraphs.begin(), EventsFor2DGraphs.end(), evId);
	if (findPosition != EventsFor2DGraphs.end()) {
		make2DGraphs = true;
		EventsFor2DGraphs.erase(findPosition);
	}

	//initialize new fit counters in case this is a new run:
	if (successfulFitCounter.find(run) == successfulFitCounter.end()) 
		successfulFitCounter[run] = failedFitCounter[run] = 0;

	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);

	//opening Clusters (made from all, closest 7, closest 9)
	edm::Handle<reco::HGCalTBClusterCollection> clusters;
  edm::Handle<reco::HGCalTBClusterCollection> clusters7;
  edm::Handle<reco::HGCalTBClusterCollection> clusters19;
	event.getByToken(HGCalTBClusterCollection_Token, clusters);
	event.getByToken(HGCalTBClusterCollection7_Token, clusters7);
	event.getByToken(HGCalTBClusterCollection19_Token, clusters19);

	//step 1: Reduce the information to energy deposits/hits in x,y per sensor/layer 
	//fill the rechits:
	for(auto Rechit : *Rechits) {	
		layer = (Rechit.id()).layer();
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap();
			Sensors[layer]->setPedestalThreshold(pedestalThreshold);
			Sensors[layer]->setZ(Layer_Z_Positions[layer], Layer_Z_X0s[layer]);	//first argument: real positon in cm, second argument: position in radiation lengths
			Sensors[layer]->setADCPerMIP(ADC_per_MIP[layer-1]);
			Sensors[layer]->setSensorSize(SensorSize);
		}
		Sensors[layer]->addHit(Rechit);
	}

	//fill the hits from the cluster collections 
	for( auto cluster : *clusters ){
		layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->registerClusterHit((*it).first, -1);
		}
	}
	for( auto cluster : *clusters7 ){
		layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->registerClusterHit((*it).first, 7);
		}
	}
	for( auto cluster : *clusters19 ){
		layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->registerClusterHit((*it).first, 19);
		}
	}

	//Event selection: sum of energies of all cells(=hits) from RecHits Collection and Clusters only must 
	//be larger than an externally configured parameter 'totalEnergyThreshold' (in MIP)
	sumEnergy = 0., sumClusterEnergy = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {	
		sumEnergy += it->second->getTotalEnergy();
		sumClusterEnergy += it->second->getTotalClusterEnergy(-1);
	}
	
	if(sumEnergy < totalEnergyThreshold) HitsVetoCounter+=1;	//make this energy dependent!
	if(sumClusterEnergy < totalEnergyThreshold) ClusterVetoCounter+=1;
	if(sumEnergy < totalEnergyThreshold && sumClusterEnergy < totalEnergyThreshold) {
		CommonVetoCounter+=1;
		return;
	} 



	//step 2: calculate impact point with technique indicated as the argument
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//subtract common first
		it->second->subtractCM();

		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(considerationMethod, weightingMethod);
		
	}


	//step 3: fill particle tracks
	std::map<int, ParticleTrack*> Tracks; 	//the integer index indicates which layer is omitted in the track calculation
	for (int i=1; i<=nLayers; i++) {
		Tracks[i] = new ParticleTrack();
		for (int j=1; j<=nLayers; j++) {
			if (i==j) continue;
			Tracks[i]->addFitPoint(Sensors[j]);
		}
		Tracks[i]->weightFitPoints(fitPointWeightingMethod);
		Tracks[i]->fitTrack(fittingMethod);
	}
	

	//step 4: calculate the deviations between each fit missing one layer and exactly that layer's true central position
	for (layer=1; layer<=nLayers; layer++) {
		layerZ_cm = Sensors[layer]->getZ_cm();
		layerZ_X0 = Sensors[layer]->getZ_X0();

		std::pair<double, double> position_predicted = Tracks[layer]->calculatePositionXY(layerZ_cm);
		x_predicted = position_predicted.first;
		y_predicted = position_predicted.second;

		std::pair<double, double> position_predicted_to_closest_cell = Sensors[layer]->getCenterOfClosestCell(position_predicted);
		x_predicted_to_closest_cell = position_predicted_to_closest_cell.first;
		y_predicted_to_closest_cell = position_predicted_to_closest_cell.second;

		std::pair<double, double> position_error_predicted = Tracks[layer]->calculatePositionErrorXY(layerZ_cm);
		x_predicted_err = position_error_predicted.first;
		y_predicted_err = position_error_predicted.second;

		if (!(x_predicted!=0 || y_predicted!=0 || x_predicted_err!=0 || y_predicted_err!=0))	{
			//default fitting has been applied, i.e. the regular fit has failed or the selected method is not implemented
			failedFitCounter[run]++;
			continue; 	//ignore those cases but count them
		}
		successfulFitCounter[run]++; 
		
		std::pair<double, double> position_true = Sensors[layer]->getCenterPosition();
		x_true = position_true.first;
		y_true = position_true.second;

		std::pair<double, double> position_true_to_closest_cell = Sensors[layer]->getCenterOfClosestCell(position_true);
		x_true_to_closest_cell = position_true_to_closest_cell.first;
		y_true_to_closest_cell = position_true_to_closest_cell.second;

		std::pair<double, double> position_error_true = Sensors[layer]->getCenterPositionError();
		x_true_err = position_error_true.first;
		y_true_err = position_error_true.second;

		deltaX = x_predicted - x_true;
		deltaY = y_predicted - y_true;
		deviation  = sqrt( pow(deltaX, 2) + pow(deltaY, 2) );

		sumFitWeights = Tracks[layer]->getSumOfEnergies();
		layerWeight = Sensors[layer]->getTotalWeight();

		//DEBUG
		if (deviation > 1000.) 
			std::cout<<"   layer: "<<layer<<"   x:  "<<x_predicted<<" - "<<x_true<<"     "<<y_predicted<<" - "<<y_true<<std::endl;
		//END OF DEBUG
		
		//fill the tree
		outTree->Fill();
		
		//update the deviations maxima on the fly
		if (make2DGraphs) {
			//store for the two 2D graphs that are written per event
			x_predicted_v.push_back(x_predicted);
			y_predicted_v.push_back(y_predicted);
			x_true_v.push_back(x_true);
			y_true_v.push_back(y_true);
			layerZ_cm_v.push_back(layerZ_cm);
		}
	}

	if (make2DGraphs) {
		std::string graphIdentifier = "run_" + std::to_string(run) + "event_" + std::to_string(evId);
		fs->make<TGraph2D>(("predicted_points_" + graphIdentifier).c_str(), "", layerZ_cm_v.size(), &(x_predicted_v[0]), &(y_predicted_v[0]), &(layerZ_cm_v[0]));
		fs->make<TGraph2D>(("true_points_" + graphIdentifier).c_str(), "", layerZ_cm_v.size(), &(x_true_v[0]), &(y_true_v[0]), &(layerZ_cm_v[0]));
		x_predicted_v.clear(); y_predicted_v.clear(); x_true_v.clear(); y_true_v.clear(); layerZ_cm_v.clear();
	}

	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	
	for(std::map<int, ParticleTrack*>::iterator it=Tracks.begin(); it!=Tracks.end(); it++) {
		delete (*it).second;
	}; Tracks.clear();
	

}// analyze ends here

void Position_Resolution_Analyzer::beginJob() {	
}

void Position_Resolution_Analyzer::endJob() {
	std::cout<<"ClusterVetos: "<<ClusterVetoCounter<<std::endl;
	std::cout<<"HitsVetos: "<<HitsVetoCounter<<std::endl;
	std::cout<<"CommonVetos: "<<CommonVetoCounter<<std::endl;
	
}

void Position_Resolution_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Position_Resolution_Analyzer);