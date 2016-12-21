/* 
 * Write out residuals and derivatives to pass forward to millepede.
 */

/**
	@Author: Thorben Quast <tquast>
		21 Dec 2016
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

#include "TFile.h"
#include "TTree.h"
  
//configuration1:
double _config1Positions[] = {0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92}; 	 //z-coordinate in cm

//configuration2:
double _config2Positions[] = {0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4}; 				//z-coordinate in cm
                     
class MillepedeNTupelizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit MillepedeNTupelizer(const edm::ParameterSet&);
		~MillepedeNTupelizer();
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

		double pedestalThreshold;
		std::vector<double> Layer_Z_Positions;
		std::vector<double> ADC_per_MIP;
		int LayersConfig;
		int nLayers;
		int SensorSize;

		double totalEnergyThreshold;

		int ClusterVetoCounter;
		int HitsVetoCounter;
		int CommonVetoCounter;

		//helper variables that are set within the event loop, i.e. are defined per event
		std::map<int, SensorHitMap*> Sensors;
		ParticleTrack* Track;

		//stuff to be written to the tree
		TTree* outTree;
		int configuration, evId, eventCounter, run; 	//eventCounter: counts the events in this analysis run to match information within ove event to each other
		double energy;
		
		//initial global parameters, corresponding to obtained value from preceding iteration
		//placeholders, should be defined per layer in the according class!
		double d_alpha, d_beta, d_gamma, d_x0, d_y0, d_z0; 

		double** residuals;
		double** errors;
		double*** derivatives;
};

MillepedeNTupelizer::MillepedeNTupelizer(const edm::ParameterSet& iConfig) {
	
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
		Layer_Z_Positions = std::vector<double>(_config1Positions, _config1Positions + sizeof(_config1Positions)/sizeof(double));
	} if (LayersConfig == 2) {
		Layer_Z_Positions = std::vector<double>(_config2Positions, _config2Positions + sizeof(_config2Positions)/sizeof(double));
	} else {
		Layer_Z_Positions = std::vector<double>(_config1Positions, _config1Positions + sizeof(_config1Positions)/sizeof(double));
	}

	pedestalThreshold = iConfig.getParameter<double>("pedestalThreshold");
	SensorSize = iConfig.getParameter<int>("SensorSize");
	nLayers = iConfig.getParameter<int>("nLayers");
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");

	totalEnergyThreshold = iConfig.getParameter<double>("totalEnergyThreshold");

	eventCounter = 0;
	residuals = new double*[nLayers];
	errors = new double*[nLayers];
	derivatives = new double**[nLayers];
	for (int i=0; i<nLayers; i++) {
		residuals[i] = new double[3];
		errors[i] = new double[3];
		derivatives[i] = new double*[3];
		for (int j=0; j<3; j++)
			derivatives[i][j] = new double[10];
	}

	//initialize tree and set Branch addresses
	outTree = fs->make<TTree>("millepede", "millepede");
	outTree->Branch("configuration", &configuration, "configuration/I");
	outTree->Branch("eventId", &evId, "eventId/I");	//event ID as it comes from the reader
	outTree->Branch("eventCounter", &eventCounter, "eventCounter/I");	//event counter, current iteration of this analysis w.r.t. the individual events
	outTree->Branch("run", &run, "run/I");
	outTree->Branch("energy", &energy, "energy/D");	//electron energy in GeV s
	
	std::string coordinates[3] = {"x", "y", "z"};
	//1-4 are local parameters
	std::string parameters[10] = {"mx", "bx", "my", "by", "dalpha", "dbeta", "dgamma", "dx0", "dy0", "dz0"};
	for (int i=0; i<nLayers; i++) {
		for (int j=0; j<3; j++){
			outTree->Branch(("residuals_layer_"+std::to_string(i+1)+"_"+coordinates[j]).c_str(), &residuals[i][j]);
			outTree->Branch(("errors_layer_"+std::to_string(i+1)+"_"+coordinates[j]).c_str(), &errors[i][j]);
			for (int k=0; k<10; k++) {
				outTree->Branch(("derivatives_layer_"+std::to_string(i+1)+"_"+coordinates[j]+"_"+parameters[k]).c_str(), &derivatives[i][j][k]);
			}
		}
	}


	ClusterVetoCounter = 0;
	HitsVetoCounter = 0;
	CommonVetoCounter = 0;


	d_alpha = 0., d_beta = 0., d_gamma = 0., d_x0 = 0., d_y0 = 0., d_z0 = 0.; 
	//possibly: read from millepede output file if procedure is iterated

}//constructor ends here

MillepedeNTupelizer::~MillepedeNTupelizer() {
	return;
}

// ------------ method called for each event  ------------
void MillepedeNTupelizer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
	edm::Handle<RunData> rd;

 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	configuration = rd->configuration;
	evId = event.id().event();
	eventCounter++;
	run = rd->run;
	energy = rd->energy;
	if (run == -1) {
		std::cout<<"Run is not in configuration file - is ignored."<<std::endl;
		return;
	}

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
		int layer = (Rechit.id()).layer();
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap();
			Sensors[layer]->setPedestalThreshold(pedestalThreshold);
			Sensors[layer]->setLabZ(Layer_Z_Positions[layer], 0.);	//first argument: real positon as measured (not aligned) in cm, second argument: position in radiation lengths
			Sensors[layer]->setAlignmentParameters(0., 0., 0., 0., 0., 0.);			//Todo: read from file or similar
			Sensors[layer]->setADCPerMIP(ADC_per_MIP[layer-1]);
			Sensors[layer]->setSensorSize(SensorSize);
		}
		Sensors[layer]->addHit(Rechit);
	}

	//fill the hits from the cluster collections 
	for( auto cluster : *clusters ){
		int layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->registerClusterHit((*it).first, -1);
		}
	}
	for( auto cluster : *clusters7 ){
		int layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->registerClusterHit((*it).first, 7);
		}
	}
	for( auto cluster : *clusters19 ){
		int layer = cluster.layer();
		for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
			Sensors[layer]->registerClusterHit((*it).first, 19);
		}
	}

	//Event selection: sum of energies of all cells(=hits) from RecHits Collection and Clusters only must 
	//be larger than an externally configured parameter 'totalEnergyThreshold' (in MIP)
	double sumEnergy = 0., sumClusterEnergy = 0.;
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
	Track = new ParticleTrack();
	for (int i=1; i<=nLayers; i++) {
		Track->addFitPoint(Sensors[i]);
	}
	Track->weightFitPoints(fitPointWeightingMethod);
	Track->fitTrack(fittingMethod);
	

	//step 4: calculate the deviations between each fit missing one layer and exactly that layer's true central position
	for (int layer=1; layer<=nLayers; layer++) {
		double layer_labZ = Sensors[layer]->getLabZ();
		double intrinsic_z = Sensors[layer]->getIntrinsicHitZPosition();	
		
		std::pair<double, double> position_predicted = Track->calculatePositionXY(layer_labZ+intrinsic_z);
		double x_predicted = position_predicted.first;
		double y_predicted = position_predicted.second;

		std::pair<double, double> position_true = Sensors[layer]->getLabHitPosition();	//should be distinguished into intrinsic position (in frame of sensor) and the lab frame
																																									  //later will require to perform according transformation with angles and displacements 
		double x_true = position_true.first;
		double y_true = position_true.second;
															
															
		residuals[layer-1][0] = x_true - x_predicted;
		residuals[layer-1][1] = y_true - y_predicted;
		residuals[layer-1][2] = 0.0;

		errors[layer-1][0] = 1.2/sqrt(12.);
		errors[layer-1][1] = 1.2/sqrt(12.);
		errors[layer-1][2] = 0.1;
		
		derivatives[layer-1][0][0] = intrinsic_z + layer_labZ;
		derivatives[layer-1][0][1] = 1.;
		derivatives[layer-1][0][2] = (intrinsic_z + layer_labZ)*d_alpha;
		derivatives[layer-1][0][3] = d_alpha;
		derivatives[layer-1][0][4] = y_predicted;
		derivatives[layer-1][0][5] = 0.;
		derivatives[layer-1][0][6] = intrinsic_z;
		derivatives[layer-1][0][7] = 1.;
		derivatives[layer-1][0][8] = 0.;
		derivatives[layer-1][0][9] = 0.;

		derivatives[layer-1][1][0] = -d_alpha*(intrinsic_z + layer_labZ);
		derivatives[layer-1][1][1] = -d_alpha;
		derivatives[layer-1][1][2] = intrinsic_z + layer_labZ;
		derivatives[layer-1][1][3] = 1.;
		derivatives[layer-1][1][4] = -x_predicted;
		derivatives[layer-1][1][5] = intrinsic_z;
		derivatives[layer-1][1][6] = 0.;
		derivatives[layer-1][1][7] = 0.;
		derivatives[layer-1][1][8] = 1.;
		derivatives[layer-1][1][9] = 0.;	

		derivatives[layer-1][2][0] = -d_gamma*(intrinsic_z + layer_labZ);
		derivatives[layer-1][2][1] = -d_gamma;
		derivatives[layer-1][2][2] = -d_beta*(intrinsic_z + layer_labZ);
		derivatives[layer-1][2][3] = -d_beta;
		derivatives[layer-1][2][4] = 0.;
		derivatives[layer-1][2][5] = -y_predicted;
		derivatives[layer-1][2][6] = -x_predicted;
		derivatives[layer-1][2][7] = 0.;
		derivatives[layer-1][2][8] = 0.;
		derivatives[layer-1][2][9] = 1.;	
		
	}
	outTree->Fill();

	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	delete Track;

}// analyze ends here

void MillepedeNTupelizer::beginJob() {	
}

void MillepedeNTupelizer::endJob() {
	std::cout<<"ClusterVetos: "<<ClusterVetoCounter<<std::endl;
	std::cout<<"HitsVetos: "<<HitsVetoCounter<<std::endl;
	std::cout<<"CommonVetos: "<<CommonVetoCounter<<std::endl;
	
}

void MillepedeNTupelizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MillepedeNTupelizer);