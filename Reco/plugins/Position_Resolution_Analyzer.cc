/* 
 * Determination of the position resolution of the setup.
 */

/**
	@Author: Thorben Quast <tquast>
		20 Febr 2017
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
#include "HGCal/DataFormats/interface/HGCalTBMultiWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Tracks.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "TFile.h"
#include "TTree.h"
  
//configuration1:
double config1Positions[] = {0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92};    //z-coordinate in cm, 1cm added to consider absorber in front of first sensor    
double config1X0Depths[] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42}; //in radiation lengths, copied from layerSumAnalyzer
//configuration2:
double config2Positions[] = {0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4};         //z-coordinate in cm, 1cm added to consider absorber in front of first sensor    
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

		struct {
  		HGCalElectronicsMap emap_;
		} essource_;
		std::string _e_mapFile;

		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;
	 	edm::EDGetToken HGCalTBClusterCollection_Token;
  		edm::EDGetToken HGCalTBClusterCollection7_Token;
  		edm::EDGetToken HGCalTBClusterCollection19_Token;

		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<MultiWireChambers> MWCToken;
		
		AlignmentParameters* alignmentParameters; //all entries are set to zero if no valid file is given 

		ConsiderationMethod considerationMethod;
		WeightingMethod weightingMethod;
		TrackFittingMethod fittingMethod;		
		FitPointWeightingMethod fitPointWeightingMethod;

		double pedestalThreshold;
		std::vector<double> Layer_Z_Positions;
		std::vector<double> Layer_Z_X0s;
		std::vector<double> ADC_per_MIP;
		int LayersConfig;
		int SensorSize;
		int nLayers;

		bool useMWCReference;


		int ClusterVetoCounter;
		int HitsVetoCounter;
		int CommonVetoCounter;

		std::map<int, int> successfulFitCounter, failedFitCounter;

		//helper variables that are set within the event loop, i.e. are defined per event
		std::map<int, SensorHitMap*> Sensors;
		std::map<int, ParticleTrack*> Tracks;

		//stuff to be written to the tree
		TTree* outTree;
		int configuration, evId, eventCounter, run, layer; 	//eventCounter: counts the events in this analysis run to match information within ove event to each other
		double energy;
		double layerWeight, layerEnergy, layerClusterEnergy, sumFitWeights, sumEnergy, sumClusterEnergy, CM_cells_count, CM_sum;
		double x_predicted, x_predicted_err, y_predicted, y_predicted_err, x_true, x_true_err, y_true, y_true_err, deltaX, deltaY;
		double x_predicted_to_closest_cell, y_predicted_to_closest_cell, x_true_to_closest_cell, y_true_to_closest_cell, layerZ_cm, layerZ_X0, deviation;

		//averaged information up to the corresponding layers
		double average_x_predicted, average_y_predicted, average_x_true, average_y_true, average_deltaX, average_deltaY; 

		//MWCs
		int useMWC;
		double MWC_x1, MWC_y1, MWC_z1, MWC_x2, MWC_y2, MWC_z2;

		std::pair<int, double> CM_tmp;	//will write the subtract_CM() return values for each layer
};

Position_Resolution_Analyzer::Position_Resolution_Analyzer(const edm::ParameterSet& iConfig) {	
	
	// initialization
	_e_mapFile = iConfig.getParameter<std::string>("e_mapFile_CERN");	
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_e_mapFile);
 	
  if (!io.load(fip.fullPath(), essource_.emap_)) {
	  throw cms::Exception("Unable to load electronics map");
	};

	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	MWCToken= consumes<MultiWireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
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
	else if (methodString == "logWeighting_2.0_1.0")
		weightingMethod = LOGWEIGHTING_20_10;
	else if (methodString == "logWeighting_2.1_1.0")
		weightingMethod = LOGWEIGHTING_21_10;
	else if (methodString == "logWeighting_2.2_1.0")
		weightingMethod = LOGWEIGHTING_22_10;
	else if (methodString == "logWeighting_2.3_1.0")
		weightingMethod = LOGWEIGHTING_23_10;
	else if (methodString == "logWeighting_2.4_1.0")
		weightingMethod = LOGWEIGHTING_24_10;
	else if (methodString == "logWeighting_2.5_1.0")
		weightingMethod = LOGWEIGHTING_25_10;
	else if (methodString == "logWeighting_2.6_1.0")
		weightingMethod = LOGWEIGHTING_26_10;
	else if (methodString == "logWeighting_2.7_1.0")
		weightingMethod = LOGWEIGHTING_27_10;
	else if (methodString == "logWeighting_2.8_1.0")
		weightingMethod = LOGWEIGHTING_28_10;
	else if (methodString == "logWeighting_2.9_1.0")
		weightingMethod = LOGWEIGHTING_29_10;
	else if (methodString == "logWeighting_3.0_1.0")
		weightingMethod = LOGWEIGHTING_30_10;
	else if (methodString == "logWeighting_3.1_1.0")
		weightingMethod = LOGWEIGHTING_31_10;
	else if (methodString == "logWeighting_3.2_1.0")
		weightingMethod = LOGWEIGHTING_32_10;
	else if (methodString == "logWeighting_3.3_1.0")
		weightingMethod = LOGWEIGHTING_33_10;
	else if (methodString == "logWeighting_3.4_1.0")
		weightingMethod = LOGWEIGHTING_34_10;
	else if (methodString == "logWeighting_3.5_1.0")
		weightingMethod = LOGWEIGHTING_35_10;
	else if (methodString == "logWeighting_3.6_1.0")
		weightingMethod = LOGWEIGHTING_36_10;
	else if (methodString == "logWeighting_3.7_1.0")
		weightingMethod = LOGWEIGHTING_37_10;
	else if (methodString == "logWeighting_3.8_1.0")
		weightingMethod = LOGWEIGHTING_38_10;
	else if (methodString == "logWeighting_3.9_1.0")
		weightingMethod = LOGWEIGHTING_39_10;
	else if (methodString == "logWeighting_4.0_1.0")
		weightingMethod = LOGWEIGHTING_40_10;
	else if (methodString == "logWeighting_4.1_1.0")
		weightingMethod = LOGWEIGHTING_41_10;
	else if (methodString == "logWeighting_4.2_1.0")
		weightingMethod = LOGWEIGHTING_42_10;
	else if (methodString == "logWeighting_4.3_1.0")
		weightingMethod = LOGWEIGHTING_43_10;
	else if (methodString == "logWeighting_4.4_1.0")
		weightingMethod = LOGWEIGHTING_44_10;
	else if (methodString == "logWeighting_4.5_1.0")
		weightingMethod = LOGWEIGHTING_45_10;
	else if (methodString == "logWeighting_4.6_1.0")
		weightingMethod = LOGWEIGHTING_46_10;
	else if (methodString == "logWeighting_4.7_1.0")
		weightingMethod = LOGWEIGHTING_47_10;
	else if (methodString == "logWeighting_4.8_1.0")
		weightingMethod = LOGWEIGHTING_48_10;
	else if (methodString == "logWeighting_4.9_1.0")
		weightingMethod = LOGWEIGHTING_49_10;
	else if (methodString == "logWeighting_5.0_1.0")
		weightingMethod = LOGWEIGHTING_50_10;
	else if (methodString == "logWeighting_6.0_1.0")
		weightingMethod = LOGWEIGHTING_60_10;
	else if (methodString == "logWeighting_7.0_1.0")
		weightingMethod = LOGWEIGHTING_70_10;
	else if (methodString == "logWeighting_2.05_1.0")
		weightingMethod = LOGWEIGHTING_205_10;
	else if (methodString == "logWeighting_2.15_1.0")
		weightingMethod = LOGWEIGHTING_215_10;
	else if (methodString == "logWeighting_2.25_1.0")
		weightingMethod = LOGWEIGHTING_225_10;
	else if (methodString == "logWeighting_2.35_1.0")
		weightingMethod = LOGWEIGHTING_235_10;
	else if (methodString == "logWeighting_2.45_1.0")
		weightingMethod = LOGWEIGHTING_245_10;
	else if (methodString == "logWeighting_2.55_1.0")
		weightingMethod = LOGWEIGHTING_255_10;
	else if (methodString == "logWeighting_2.65_1.0")
		weightingMethod = LOGWEIGHTING_265_10;
	else if (methodString == "logWeighting_2.75_1.0")
		weightingMethod = LOGWEIGHTING_275_10;
	else if (methodString == "logWeighting_2.85_1.0")
		weightingMethod = LOGWEIGHTING_285_10;
	else if (methodString == "logWeighting_2.95_1.0")
		weightingMethod = LOGWEIGHTING_295_10;
	else if (methodString == "logWeighting_3.05_1.0")
		weightingMethod = LOGWEIGHTING_305_10;
	else if (methodString == "logWeighting_3.15_1.0")
		weightingMethod = LOGWEIGHTING_315_10;
	else if (methodString == "logWeighting_3.25_1.0")
		weightingMethod = LOGWEIGHTING_325_10;
	else if (methodString == "logWeighting_3.35_1.0")
		weightingMethod = LOGWEIGHTING_335_10;
	else if (methodString == "logWeighting_3.45_1.0")
		weightingMethod = LOGWEIGHTING_345_10;
	else if (methodString == "logWeighting_3.55_1.0")
		weightingMethod = LOGWEIGHTING_355_10;
	else if (methodString == "logWeighting_3.65_1.0")
		weightingMethod = LOGWEIGHTING_365_10;
	else if (methodString == "logWeighting_3.75_1.0")
		weightingMethod = LOGWEIGHTING_375_10;
	else if (methodString == "logWeighting_3.85_1.0")
		weightingMethod = LOGWEIGHTING_385_10;
	else if (methodString == "logWeighting_3.95_1.0")
		weightingMethod = LOGWEIGHTING_395_10;
	else if (methodString == "logWeighting_4.05_1.0")
		weightingMethod = LOGWEIGHTING_405_10;
	else if (methodString == "logWeighting_4.15_1.0")
		weightingMethod = LOGWEIGHTING_415_10;
	else if (methodString == "logWeighting_4.25_1.0")
		weightingMethod = LOGWEIGHTING_425_10;
	else if (methodString == "logWeighting_4.35_1.0")
		weightingMethod = LOGWEIGHTING_435_10;
	else if (methodString == "logWeighting_4.45_1.0")
		weightingMethod = LOGWEIGHTING_445_10;
	else if (methodString == "logWeighting_4.55_1.0")
		weightingMethod = LOGWEIGHTING_455_10;
	else if (methodString == "logWeighting_4.65_1.0")
		weightingMethod = LOGWEIGHTING_465_10;
	else if (methodString == "logWeighting_4.75_1.0")
		weightingMethod = LOGWEIGHTING_475_10;
	else if (methodString == "logWeighting_4.85_1.0")
		weightingMethod = LOGWEIGHTING_485_10;
	else if (methodString == "logWeighting_4.95_1.0")
		weightingMethod = LOGWEIGHTING_495_10;
	else 
		weightingMethod = DEFAULTWEIGHTING;

	//read the track fitting method
	methodString = iConfig.getParameter<std::string>("fittingMethod");
	if (methodString == "lineAnalytical")
		fittingMethod = LINEFITANALYTICAL;
	else if (methodString == "lineTGraphErrors")
		fittingMethod = LINEFITTGRAPHERRORS;
	else if (methodString == "pol2TGraphErrors")
		fittingMethod = POL2TGRAPHERRORS;
	else if (methodString == "pol3TGraphErrors")
		fittingMethod = POL3TGRAPHERRORS;
	else if (methodString == "gblTrack")
		fittingMethod = GBLTRACK;
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

	eventCounter = 0;

	pedestalThreshold = iConfig.getParameter<double>("pedestalThreshold");
	SensorSize = iConfig.getParameter<int>("SensorSize");
	nLayers = iConfig.getParameter<int>("nLayers");
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");


	useMWCReference = iConfig.getParameter<bool>("useMWCReference");

	//initialize tree and set Branch addresses
	outTree = fs->make<TTree>("deviations", "deviations");
	outTree->Branch("configuration", &configuration, "configuration/I");
	outTree->Branch("eventId", &evId, "eventId/I");	//event ID as it comes from the reader, as it is stored in the txt files
	outTree->Branch("eventCounter", &eventCounter, "eventCounter/I");	//event counter, current iteration, indexing occurs chronologically in the readin plugins
	outTree->Branch("run", &run, "run/I");
	outTree->Branch("layer", &layer, "layer/I");
	outTree->Branch("energy", &energy, "energy/D");	//electron energy in GeV
	
	outTree->Branch("CM_sum", &CM_sum, "CM_sum/D");
	outTree->Branch("CM_cells_count", &CM_cells_count, "CM_cells_count/D");
	outTree->Branch("layerEnergy", &layerEnergy, "layerEnergy/D");
	outTree->Branch("layerClusterEnergy", &layerClusterEnergy, "layerClusterEnergy/D");
	outTree->Branch("layerWeight", &layerWeight, "layerWeight/D");
	outTree->Branch("sumEnergy", &sumEnergy, "sumEnergy/D");
	outTree->Branch("sumClusterEnergy", &sumClusterEnergy, "sumClusterEnergy/D");
	outTree->Branch("sumFitWeights", &sumFitWeights, "sumFitWeights/D");
	
	outTree->Branch("x_true", &x_true, "x_true/D");
	outTree->Branch("x_true_to_closest_cell", &x_true_to_closest_cell, "x_true_to_closest_cell/D");
	outTree->Branch("x_true_err", &x_true_err, "x_true_err/D");
	outTree->Branch("y_true", &y_true, "y_true/D");
	outTree->Branch("y_true_to_closest_cell", &y_true_to_closest_cell, "y_true_to_closest_cell/D");
	outTree->Branch("y_true_err", &y_true_err, "y_true_err/D");
	
	outTree->Branch("x_predicted", &x_predicted, "x_predicted/D");
	outTree->Branch("x_predicted_to_closest_cell", &x_predicted_to_closest_cell, "x_predicted_to_closest_cell/D");
	outTree->Branch("x_predicted_err", &x_predicted_err, "x_predicted_err/D");
	outTree->Branch("y_predicted", &y_predicted, "y_predicted/D");
	outTree->Branch("y_predicted_to_closest_cell", &y_predicted_to_closest_cell, "y_predicted_to_closest_cell/D");
	outTree->Branch("y_predicted_err", &y_predicted_err, "y_predicted_err/D");
	
	outTree->Branch("useMWC", &useMWC, "uswMWC/I");
	outTree->Branch("MWC_x1", &MWC_x1, "MWC_x1/D");
	outTree->Branch("MWC_y1", &MWC_y1, "MWC_y1/D");
	outTree->Branch("MWC_z1", &MWC_z1, "MWC_z1/D");
	outTree->Branch("MWC_x2", &MWC_x2, "MWC_x2/D");
	outTree->Branch("MWC_y2", &MWC_y2, "MWC_y2/D");
	outTree->Branch("MWC_z2", &MWC_z2, "MWC_z2/D");
	
	outTree->Branch("deltaX", &deltaX, "deltaX/D");
	outTree->Branch("deltaY", &deltaY, "deltaY/D");
	outTree->Branch("layerZ_cm", &layerZ_cm, "layerZ_cm/D");
	outTree->Branch("layerZ_X0", &layerZ_X0, "layerZ_X0/D");
	outTree->Branch("deviation", &deviation, "deviation/D");

	outTree->Branch("average_x_predicted", &average_x_predicted, "average_x_predicted/D");
	outTree->Branch("average_y_predicted", &average_y_predicted, "average_y_predicted/D");
	outTree->Branch("average_x_true", &average_x_true, "average_x_true/D");
	outTree->Branch("average_y_true", &average_y_true, "average_y_true/D");
	outTree->Branch("average_deltaX", &average_deltaX, "average_deltaX/D");
	outTree->Branch("average_deltaY", &average_deltaY, "average_deltaY/D");

	alignmentParameters = new AlignmentParameters(iConfig.getParameter<std::vector<std::string> >("alignmentParameterFiles")); 

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
	eventCounter = rd->event;
	energy = rd->energy;
	if (rd->hasDanger) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because it has DANGER=true"<<std::endl;
		return;
	}
	if (useMWCReference && ! rd->hasValidMWCMeasurement) {
		//std::cout<<"Event "<<event.id().event()<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because it has an invalid MWC measurement"<<std::endl;
		return;	
	}
	if (run == -1) {
		std::cout<<"Run is not in configuration file - is ignored."<<std::endl;
		return;
	}


	edm::Handle<MultiWireChambers> mwcs;
	event.getByToken(MWCToken, mwcs);

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
			Sensors[layer] = new SensorHitMap(layer);
			Sensors[layer]->setPedestalThreshold(pedestalThreshold);
			Sensors[layer]->setLabZ(Layer_Z_Positions[layer-1], Layer_Z_X0s[layer-1]);	//first argument: real positon as measured (not aligned) in cm, second argument: position in radiation lengths

			Sensors[layer]->setAlignmentParameters(alignmentParameters->getValue(run, 100*layer + 21), 0.0, 0.0,
				alignmentParameters->getValue(run, 100*layer + 11), alignmentParameters->getValue(run, 100*layer + 12), 0.0);	
			Sensors[layer]->setSensorSize(SensorSize);

			double X0sum = 0;
			for (int _x = 0; _x<(int)layer; _x++) X0sum += Layer_Z_X0s[_x];
			Sensors[layer]->setParticleEnergy(energy - gblhelpers::computeEnergyLoss(X0sum, energy));
		}
		uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
		HGCalTBElectronicsId eid(EID);	 
		int skiRocIndex = (eid.iskiroc() - 1) > 0 ? eid.iskiroc() - 1 : 0;	
		Sensors[layer]->addHit(Rechit, ADC_per_MIP[skiRocIndex]);
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

	//Possible event selection: sum of energies of all cells(=hits) from RecHits Collection and Clusters
	sumEnergy = 0., sumClusterEnergy = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {	
		sumEnergy += it->second->getTotalEnergy();
		sumClusterEnergy += it->second->getTotalClusterEnergy(-1);
	}
	

	//step 2: calculate impact point with technique indicated as the argument
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//subtract common first
		CM_tmp = it->second->subtractCM();
		CM_cells_count = CM_tmp.first;
		CM_sum = CM_tmp.second;

		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(considerationMethod, weightingMethod);
	}

	//Step 3: add MWCs to the setup if useMWCReference option is set true
	if (useMWCReference) {
		useMWC = 1;

		double mwc_resolution = 0.005; //cm
		Sensors[(nLayers+1)] = new SensorHitMap((nLayers+1));				//attention: This is specifically tailored for the 8-layer setup
		Sensors[(nLayers+1)]->setLabZ(mwcs->at(0).z, 0.001);
		Sensors[(nLayers+1)]->setCenterHitPosition(mwcs->at(0).x, mwcs->at(0).y ,mwc_resolution , mwc_resolution);
		Sensors[(nLayers+1)]->setParticleEnergy(energy);
		Sensors[(nLayers+1)]->setAlignmentParameters(alignmentParameters->getValue(energy, 100*(nLayers+1) + 21), 0.0, 0.0,
				alignmentParameters->getValue(energy, 100*(nLayers+1) + 11), alignmentParameters->getValue(energy, 100*(nLayers+1) + 12), 0.0);	
		Sensors[(nLayers+1)]->setResidualResolution(mwc_resolution);	
		MWC_x1 = Sensors[nLayers+1]->getLabHitPosition().first; //mwcs->at(0).x;
		MWC_y1 = Sensors[nLayers+1]->getLabHitPosition().second; //mwcs->at(0).y;
		MWC_z1 = Sensors[nLayers+1]->getLabZ() + Sensors[nLayers+1]->getIntrinsicHitZPosition(); //mwcs->at(0).z;

		Sensors[(nLayers+2)] = new SensorHitMap((nLayers+2));				//attention: This is specifically tailored for the 8-layer setup
		Sensors[(nLayers+2)]->setLabZ(mwcs->at(1).z, 0.001);
		Sensors[(nLayers+2)]->setCenterHitPosition(mwcs->at(1).x, mwcs->at(1).y ,mwc_resolution , mwc_resolution);
		Sensors[(nLayers+2)]->setParticleEnergy(energy);
		Sensors[(nLayers+2)]->setAlignmentParameters(alignmentParameters->getValue(energy, 100*(nLayers+2) + 21), 0.0, 0.0,
				alignmentParameters->getValue(energy, 100*(nLayers+2) + 11), alignmentParameters->getValue(energy, 100*(nLayers+2) + 12), 0.0);	
		Sensors[(nLayers+2)]->setResidualResolution(mwc_resolution);
		MWC_x2 = Sensors[nLayers+2]->getLabHitPosition().first; //mwcs->at(1).x;
		MWC_y2 = Sensors[nLayers+2]->getLabHitPosition().second; //mwcs->at(1).y;
		MWC_z2 = Sensors[nLayers+2]->getLabZ() + Sensors[nLayers+1]->getIntrinsicHitZPosition(); //mwcs->at(1).z;

	} else {
		useMWC = 0;
		MWC_x1 = MWC_y1 = MWC_z1 = MWC_x2 = MWC_y2 = MWC_z2 = -999;
	}
	
	//step 4: fill particle tracks
	std::map<int, ParticleTrack*> Tracks; 	//the integer index indicates which layer is omitted in the track calculation
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		int i = it->first;
		
		Tracks[i] = new ParticleTrack();
		Tracks[i]->addReferenceSensor(Sensors[i]);

		for (std::map<int, SensorHitMap*>::iterator jt=Sensors.begin(); jt!=Sensors.end(); jt++) {
			int j = jt->first;
			if (i==j) continue;

			if (i<=nLayers) {
				if(useMWCReference && j>nLayers) Tracks[i]->addFitPoint(Sensors[j]);
				else if(!useMWCReference && j<=nLayers) Tracks[i]->addFitPoint(Sensors[j]);
			} else {
				if(useMWCReference && j<=nLayers) 
					Tracks[i]->addFitPoint(Sensors[j]);	
				
			}
		}

		Tracks[i]->weightFitPoints(fitPointWeightingMethod);
		Tracks[i]->fitTrack(fittingMethod);
	}


	//step 5: calculate the deviations between each fit missing one layer and exactly that layer's true central position
	double sum_x_predicted = 0, sum_y_predicted = 0, sum_x_true = 0, sum_y_true = 0, sum_energy = 0;

	layerZ_X0 = 0;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		layer = it->first;
		sumFitWeights = Tracks[layer]->getSumOfEnergies();
		layerEnergy = Sensors[layer]->getTotalEnergy();
		sum_energy += layerEnergy;
		layerClusterEnergy = Sensors[layer]->getTotalClusterEnergy(-1);
		layerWeight = Sensors[layer]->getTotalWeight();
		layerZ_cm = Sensors[layer]->getLabZ() + Sensors[layer]->getIntrinsicHitZPosition();
		layerZ_X0 += Sensors[layer]->getX0();
		
		std::pair<double, double> position_predicted = Tracks[layer]->calculateReferenceXY();
		x_predicted = position_predicted.first;
		sum_x_predicted  += x_predicted*layerEnergy;
		y_predicted = position_predicted.second;
		sum_y_predicted  += y_predicted*layerEnergy;

		std::pair<double, double> position_predicted_to_closest_cell = Sensors[layer]->getCenterOfClosestCell(position_predicted);
		x_predicted_to_closest_cell = position_predicted_to_closest_cell.first;
		y_predicted_to_closest_cell = position_predicted_to_closest_cell.second;
		std::pair<double, double> position_error_predicted = Tracks[layer]->calculateReferenceErrorXY();
		x_predicted_err = position_error_predicted.first;
		y_predicted_err = position_error_predicted.second;


		if (!(x_predicted!=0 || y_predicted!=0 || x_predicted_err!=0 || y_predicted_err!=0))	{
			//default fitting has been applied, i.e. the regular fit has failed or the selected method is not implemented
			failedFitCounter[run]++;
			continue; 	//ignore those cases but count them
		}
		successfulFitCounter[run]++; 

		std::pair<double, double> position_true = Sensors[layer]->getLabHitPosition();
		x_true = position_true.first;
		sum_x_true += x_true*layerEnergy;
		y_true = position_true.second;
		sum_y_true += y_true*layerEnergy;
		
		std::pair<double, double> position_true_to_closest_cell = Sensors[layer]->getCenterOfClosestCell(position_true);
		x_true_to_closest_cell = position_true_to_closest_cell.first;
		y_true_to_closest_cell = position_true_to_closest_cell.second;
		std::pair<double, double> position_error_true = Sensors[layer]->getHitPositionError();
		x_true_err = position_error_true.first;
		y_true_err = position_error_true.second;
		deltaX = x_true - x_predicted;
		deltaY = y_true - y_predicted;
		deviation  = sqrt( pow(deltaX, 2) + pow(deltaY, 2) );


		average_x_predicted = layer <= 8 ? sum_x_predicted / sum_energy : -999;
		average_y_predicted = layer <= 8 ? sum_y_predicted / sum_energy : -999;
		average_x_true = layer <= 8 ? sum_x_true / sum_energy : -999;
		average_y_true = layer <= 8 ? sum_y_true / sum_energy : -999;
		average_deltaX = layer <= 8 ? average_x_true - average_x_predicted : -999;
		average_deltaY = layer <= 8 ? average_y_true - average_y_predicted : -999;
		//DEBUG
		if (deviation > 1000.) {
			std::cout<<"Event: "<<eventCounter<<std::endl;
			std::cout<<"   layer: "<<layer<<"   x:  "<<x_predicted<<" - "<<x_true<<"     "<<y_predicted<<" - "<<y_true<<std::endl;
		}
		//END OF DEBUG
		
		//fill the tree
		outTree->Fill();
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

	delete alignmentParameters;
}

void Position_Resolution_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Position_Resolution_Analyzer);