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
  
double config1Positions[] = {0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92};
double config2Positions[] = {0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4};
                     
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

		std::vector<int> EventsFor2DGraphs;
		double pedestalThreshold;
		std::vector<double> Layer_Z_Positions;
		std::vector<double> ADC_per_MIP;
		int LayersConfig;
		int nLayers;
		int SensorSize;

		std::map<int, int> successfulFitCounter, failedFitCounter;


		//helper variables that are set within the event loop, i.e. are defined per event
		std::map<int, SensorHitMap*> Sensors;
		std::map<int, ParticleTrack*> Tracks;
		std::vector<double> x_predicted_v,y_predicted_v, x_true_v, y_true_v, layerZ_v;

		//stuff to be written to the tree
		TTree* outTree;
		int evId, run, layer;
		double energy, layerThickness;
		double x_predicted, y_predicted, x_true, y_true, deltaX, deltaY, layerZ, deviation;
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
	else if (methodString == "logWeighting_5.0_1.0")
		weightingMethod = LOGWEIGHTING_50_10;
	else if (methodString == "logWeighting_5.0_0.5")
		weightingMethod = LOGWEIGHTING_50_05;
	else if (methodString == "logWeighting_7.0_1.0")
		weightingMethod = LOGWEIGHTING_70_10;
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

	//read the layer configuration
	LayersConfig = iConfig.getParameter<int>("layers_config");
	if (LayersConfig == 1) Layer_Z_Positions = std::vector<double>(config1Positions, config1Positions + sizeof(config1Positions)/sizeof(double));
	if (LayersConfig == 2) Layer_Z_Positions = std::vector<double>(config2Positions, config2Positions + sizeof(config2Positions)/sizeof(double));
	else Layer_Z_Positions = std::vector<double>(config1Positions, config1Positions + sizeof(config1Positions)/sizeof(double));


	pedestalThreshold = iConfig.getParameter<double>("pedestalThreshold");
	SensorSize = iConfig.getParameter<int>("SensorSize");
	nLayers = iConfig.getParameter<int>("nLayers");
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");

	//making 2DGraphs per event?
	EventsFor2DGraphs = iConfig.getParameter<std::vector<int> >("EventsFor2DGraphs");

	//initialize tree and set Branch addresses
	outTree = fs->make<TTree>("deviations", "deviations");
	outTree->Branch("eventId", &evId, "eventId/I");
	outTree->Branch("run", &run, "run/I");
	outTree->Branch("layer", &layer, "layer/I");
	outTree->Branch("energy", &energy, "energy/D");
	outTree->Branch("layerThickness", &layerThickness, "layerThickness/D");
	outTree->Branch("x_predicted", &x_predicted, "x_predicted/D");
	outTree->Branch("y_predicted", &y_predicted, "y_predicted/D");
	outTree->Branch("x_true", &x_true, "x_true/D");
	outTree->Branch("y_true", &y_true, "y_true/D");
	outTree->Branch("deltaX", &deltaX, "deltaX/D");
	outTree->Branch("deltaY", &deltaY, "deltaY/D");
	outTree->Branch("layerZ", &layerZ, "layerZ/D");
	outTree->Branch("deviation", &deviation, "deviation/D");
	//	double energy, layerThickness;
	//	double x_predicted, y_predicted, x_true, y_true, deltaX, deltaY, layerZ, deviation;


}//constructor ends here

Position_Resolution_Analyzer::~Position_Resolution_Analyzer() {
	return;
}

// ------------ method called for each event  ------------
void Position_Resolution_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
	edm::Handle<RunData> rd;

 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	evId = event.id().event();
	run = rd->run;
	energy = rd->energy;
	layerThickness = rd->layerThickness;
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
			Sensors[layer]->setZ(Layer_Z_Positions[layer]);
			Sensors[layer]->setADCPerMIP(ADC_per_MIP[layer-1]);
			Sensors[layer]->setSensorSize(SensorSize);
		}
		Sensors[layer]->addHit(Rechit);
	}

	//fill the hits from the cluster collections 
	for( auto cluster : *clusters ){
		layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->addClusterHit((*it).first, -1);
		}
	}
	for( auto cluster : *clusters7 ){
		layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->addClusterHit((*it).first, 7);
		}
	}
	for( auto cluster : *clusters19 ){
		layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->addClusterHit((*it).first, 19);
		}
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
		Tracks[i]->fitTrack(fittingMethod);
	}
	
	//step 4: calculate the deviations between each fit missing one layer and exactly that layer's true central position
	for (int layer=1; layer<=nLayers; layer++) {
		layerZ = Sensors[layer]->getZ();
		x_predicted = Tracks[layer]->calculatePositionXY(layerZ).first;
		y_predicted = Tracks[layer]->calculatePositionXY(layerZ).second;
		
		if (x_predicted==0 && y_predicted==0)	{
			//default fitting has been applied, i.e. the regular fit has failed or the selected method is not implemented
			failedFitCounter[run]++;
			continue; 	//ignore those cases but count them
		}
		successfulFitCounter[run]++; 
		
		x_true = Sensors[layer]->getCenterPosition().first;
		y_true = Sensors[layer]->getCenterPosition().second;

		deltaX = x_predicted - x_true;
		deltaY = y_predicted - y_true;
		deviation  = sqrt( pow(deltaX, 2) + pow(deltaY, 2) );
		
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
			layerZ_v.push_back(layerZ);
		}
	}

	if (make2DGraphs) {
		std::string graphIdentifier = "run_" + std::to_string(run) + "event_" + std::to_string(evId);
		fs->make<TGraph2D>(("predicted_points_" + graphIdentifier).c_str(), "", layerZ_v.size(), &(x_predicted_v[0]), &(y_predicted_v[0]), &(layerZ_v[0]));
		fs->make<TGraph2D>(("true_points_" + graphIdentifier).c_str(), "", layerZ_v.size(), &(x_true_v[0]), &(y_true_v[0]), &(layerZ_v[0]));
		x_predicted_v.clear(); y_predicted_v.clear(); x_true_v.clear(); y_true_v.clear(); layerZ_v.clear();
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
	//optained information are not processed into according ROOT objects for final storage
}

void Position_Resolution_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Position_Resolution_Analyzer);