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
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Reco/interface/Mille.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Tracks.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "Alignment/ReferenceTrajectories/interface/MilleBinary.h"

/*************/
/* Some hard coded numbers:  */  

//configuration1:
double _config1Positions[] = {1.0, 6.35, 11.52, 15.44, 19.52, 20.67, 24.78, 26.92};    //z-coordinate in cm, 1cm added to consider absorber in front of first sensor    
double _config1X0Depths[] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42}; //in radiation lengths, copied from layerSumAnalyzer
//configuration2:
double _config2Positions[] = {1.0, 5.67, 10.84, 15.27, 20.25, 21.4, 26.8, 32.4};         //z-coordinate in cm, 1cm added to consider absorber in front of first sensor    
double _config2X0Depths[] = {5.048, 3.412, 3.412, 2.866, 2.512, 1.625, 2.368, 6.021}; //in radiation lengths, copied from layerSumAnalyzer

double sigma_res_x[] = {0.2, 0.15, 0.15, 0.15, 0.17, 0.18, 0.21, 0.25};					//width of residuals in x dimension (from logweighting(5,1), no reweighting of errors, all cells considered, 2MIP threshold)
double sigma_res_y[] = {0.21, 0.16, 0.16, 0.16, 0.18, 0.18, 0.21, 0.26};					//width of residuals in y dimension (from logweighting(5,1), no reweighting of errors, all cells considered, 2MIP threshold)

/*************/

class MillepedeBinaryWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit MillepedeBinaryWriter(const edm::ParameterSet&);
		~MillepedeBinaryWriter();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
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
		std::vector<double> Layer_Z_X0s;
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

		int run; 
		double energy;
		
		Mille* mille;
		gbl::MilleBinary* milleBinary;

  	int NLC, NGLperLayer, NGL;
		float rMeas, sigma;
		float *derLc, *derGl;
		int *label;
};

MillepedeBinaryWriter::MillepedeBinaryWriter(const edm::ParameterSet& iConfig) {
	
	// initialization
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
		Layer_Z_Positions = std::vector<double>(_config1Positions, _config1Positions + sizeof(_config1Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(_config1X0Depths, _config1X0Depths + sizeof(_config1X0Depths)/sizeof(double));
	} if (LayersConfig == 2) {
		Layer_Z_Positions = std::vector<double>(_config2Positions, _config2Positions + sizeof(_config2Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(_config2X0Depths, _config2X0Depths + sizeof(_config2X0Depths)/sizeof(double));
	} else {
		Layer_Z_Positions = std::vector<double>(_config1Positions, _config1Positions + sizeof(_config1Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(_config1X0Depths, _config1X0Depths + sizeof(_config1X0Depths)/sizeof(double));
	}

	pedestalThreshold = iConfig.getParameter<double>("pedestalThreshold");
	SensorSize = iConfig.getParameter<int>("SensorSize");
	nLayers = iConfig.getParameter<int>("nLayers");
	ADC_per_MIP = iConfig.getParameter<std::vector<double> >("ADC_per_MIP");

	totalEnergyThreshold = iConfig.getParameter<double>("totalEnergyThreshold");
	
	if (fittingMethod==GBLTRACK) {
		milleBinary = new gbl::MilleBinary((iConfig.getParameter<std::string>("binaryFile")).c_str());
		mille = NULL;
	} else{
		milleBinary = NULL;
		mille = new Mille((iConfig.getParameter<std::string>("binaryFile")).c_str());
	}
  

  NLC = 4;
  NGLperLayer = 3;
  NGL = nLayers*NGLperLayer;
	rMeas = 0.;
	sigma = 0.;
	derLc = new float[NLC];
	derGl = new float[NGL];
	label = new int[NGL];

	for (int l=1; l<=nLayers; l++) {
		label[(l-1)*NGLperLayer + 0] = l*100 + 11;
		label[(l-1)*NGLperLayer + 1] = l*100 + 12;
		label[(l-1)*NGLperLayer + 2] = l*100 + 21;
	}

	ClusterVetoCounter = 0;
	HitsVetoCounter = 0;
	CommonVetoCounter = 0;

}//constructor ends here

MillepedeBinaryWriter::~MillepedeBinaryWriter() {
	return;
}

// ------------ method called for each event  ------------
void MillepedeBinaryWriter::analyze(const edm::Event& event, const edm::EventSetup& setup) {
	edm::Handle<RunData> rd;

 	//get the relevant event information
	event.getByToken(RunDataToken, rd);

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
	for (auto Rechit : *Rechits) {	
		int layer = (Rechit.id()).layer();
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap(layer);
			Sensors[layer]->setPedestalThreshold(pedestalThreshold);
			Sensors[layer]->setParticleEnergy(energy);
			Sensors[layer]->setLabZ(Layer_Z_Positions[layer-1], Layer_Z_X0s[layer-1]);	//first argument: real positon as measured (not aligned) in cm, second argument: position in radiation lengths
			Sensors[layer]->setAlignmentParameters(0.0, 0.0, 0.0,
				0.0, 0.0, 0.0);	
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
	
	if (fittingMethod==GBLTRACK) {
		Track->gblTrackToMilleBinary(milleBinary);
	}else{
		//step 4: calculate the deviations between each fit missing one layer and exactly that layer's true central position
		for (int layer=1; layer<=nLayers; layer++) {
			double layer_labZ = Sensors[layer]->getLabZ();
			double intrinsic_z = Sensors[layer]->getIntrinsicHitZPosition();	
			
			std::pair<double, double> position_predicted = Track->calculatePositionXY(layer_labZ+intrinsic_z, layer);
			double x_predicted = position_predicted.first;
			double y_predicted = position_predicted.second;

			//std::pair<double, double> position_predicted_err = Track->calculatePositionErrorXY(layer_labZ+intrinsic_z);
			//double x_predicted_err = position_predicted_err.first;
			//double y_predicted_err = position_predicted_err.second;

			std::pair<double, double> position_true = Sensors[layer]->getHitPosition();	
			double x_true = position_true.first;
			double y_true = position_true.second;
			
			//std::pair<double, double> position_true_err = Sensors[layer]->getHitPositionError();	
			Sensors[layer]->getHitPositionError();	
			//double x_true_err = position_true_err.first;
			//double y_true_err = position_true_err.second;

		//step5: calculate the necessary derivatives for Mille and fill them into the binary file			
			//reset all global parameters:
			for (int k=0; k<NGL; k++){
				derGl[k] = 0.;
			}

			//the x-coordinate								
			derLc[0] = layer_labZ;
			derLc[1] = 1.;
			derLc[2] = 0.;
			derLc[3] = 0.;
			
			derGl[(layer-1)*NGLperLayer+0] = 1.;		
			derGl[(layer-1)*NGLperLayer+1] = 0.;		

			rMeas = x_true - x_predicted;
			sigma = sigma_res_x[layer-1];
			mille->mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);


			//the y-coordinate
			derLc[0] = 0.;
			derLc[1] = 0.;
			derLc[2] = layer_labZ;
			derLc[3] = 1.;

			derGl[(layer-1)*NGLperLayer+0] = 0.;		
			derGl[(layer-1)*NGLperLayer+1] = 1.;			

			rMeas = y_true - y_predicted;
			sigma = sigma_res_y[layer-1];
			mille->mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);
		}
		mille->end();
	}
	
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	delete Track;

}// analyze ends here

void MillepedeBinaryWriter::beginJob() {	
}

void MillepedeBinaryWriter::endJob() {
	if (mille != NULL) {
		mille->kill();
		delete mille;
	}
	if (milleBinary != NULL) {
		delete milleBinary;
	}
	delete derLc; delete derGl; delete label;
	std::cout<<"ClusterVetos: "<<ClusterVetoCounter<<std::endl;
	std::cout<<"HitsVetos: "<<HitsVetoCounter<<std::endl;
	std::cout<<"CommonVetos: "<<CommonVetoCounter<<std::endl;
	
}

void MillepedeBinaryWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MillepedeBinaryWriter);