/* 
 * Write out residuals and derivatives to pass forward to millepede. Relies on the fact that all layers have proper hits. 
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
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBMultiWireChamberData.h"
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
double _config1Positions[] = {0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92};    //z-coordinate in cm, 1cm added to consider absorber in front of first sensor    
double _config1X0Depths[] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42}; //in radiation lengths, copied from layerSumAnalyzer
//configuration2:
double _config2Positions[] = {0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4};         //z-coordinate in cm, 1cm added to consider absorber in front of first sensor    
double _config2X0Depths[] = {5.048, 3.412, 3.412, 2.866, 2.512, 1.625, 2.368, 6.021}; //in radiation lengths, copied from layerSumAnalyzer

double sigma_res[] = {0.205, 0.155, 0.155, 0.155, 0.175, 0.18, 0.21, 0.255};					//width of residuals in x dimension (from logweighting(5,1), no reweighting of errors, all cells considered, 2MIP threshold)

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
		struct {
  		HGCalElectronicsMap emap_;
		} essource_;
		std::string _e_mapFile;

		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;
	 	edm::EDGetToken HGCalTBClusterCollection_Token;
  		edm::EDGetToken HGCalTBClusterCollection7_Token;
  		edm::EDGetToken HGCalTBClusterCollection19_Token;

		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<MultiWireChambers> MWCToken;

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
		int nPlanes;
		int SensorSize;

		double totalEnergyThreshold;
		bool useMWCReference;

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
	_e_mapFile = iConfig.getParameter<std::string>("e_mapFile_CERN");	
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_e_mapFile);
 	
  	if (!io.load(fip.fullPath(), essource_.emap_)) {
	  throw cms::Exception("Unable to load electronics map");
	};

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

	useMWCReference = iConfig.getParameter<bool>("useMWCReference");
	
	if (fittingMethod==GBLTRACK) {
		milleBinary = new gbl::MilleBinary((iConfig.getParameter<std::string>("binaryFile")).c_str());
		mille = NULL;
	} else{
		milleBinary = NULL;
		mille = new Mille((iConfig.getParameter<std::string>("binaryFile")).c_str());
	}
  

	nPlanes = nLayers + (useMWCReference ? 2: 0);

	NLC = 4;
	NGLperLayer = 3;
	NGL = nPlanes*NGLperLayer;
	rMeas = 0.;
	sigma = 0.;
	derLc = new float[NLC];
	derGl = new float[NGL];
	label = new int[NGL];

	for (int l=1; l<=nPlanes; l++) {
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

	if (rd->hasDanger) {
		std::cout<<"Event "<<event.id().event()<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because it has DANGER=true"<<std::endl;
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

	//get the multi wire chambers
	edm::Handle<MultiWireChambers> mwcs;
	event.getByToken(MWCToken, mwcs);

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
			Sensors[layer]->setLabZ(Layer_Z_Positions[layer-1], Layer_Z_X0s[layer-1]);	//first argument: real positon as measured (not aligned) in cm, second argument: position in radiation lengths
			Sensors[layer]->setAlignmentParameters(0.0, 0.0, 0.0,
				0.0, 0.0, 0.0);	
			Sensors[layer]->setResidualResolution(sigma_res[layer-1]);
			Sensors[layer]->setSensorSize(SensorSize);
	
			uint32_t EID = essource_.emap_.detId2eid(Rechit.id());
			HGCalTBElectronicsId eid(EID);	 
			int skiRocIndex = (eid.iskiroc() - 1) > 0 ? eid.iskiroc() - 1 : 0;	
			Sensors[layer]->setADCPerMIP(ADC_per_MIP[skiRocIndex]);

			double X0sum = 0;
			for (int _x = 0; _x<(int)layer; _x++) X0sum += Layer_Z_X0s[_x];
			Sensors[layer]->setParticleEnergy(energy - gblhelpers::computeEnergyLoss(X0sum, energy));
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

	//step 3: fill all remaining layers that do not have any hits (e.g. in simulation with low energetic electrons)
	for (int layer=1; layer<=nLayers; layer++) {
		if (Sensors.find(layer)==Sensors.end()) {
			Sensors[layer] = new SensorHitMap(layer);
		
			double X0sum = 0;
			for (int _x = 0; _x<(int)layer; _x++) X0sum += Layer_Z_X0s[_x];
			Sensors[layer]->setParticleEnergy(energy - gblhelpers::computeEnergyLoss(X0sum, energy));

			Sensors[layer]->setCenterHitPosition(0, 0, 12./sqrt(12.), 12./sqrt(12.));	
		}
	}

	//Step 4: add MWCs to the setup if useMWCReference option is set true
	if (useMWCReference) {
		double mwc_resolution = 0.005; //cm
		
		Sensors[(nLayers+1)] = new SensorHitMap((nLayers+1));				//attention: This is specifically tailored for the 8-layer setup
		Sensors[(nLayers+1)]->setLabZ(mwcs->at(0).z, 0.001);
		Sensors[(nLayers+1)]->setCenterHitPosition(mwcs->at(0).x, mwcs->at(0).y ,mwc_resolution , mwc_resolution);
		Sensors[(nLayers+1)]->setParticleEnergy(energy);
		Sensors[(nLayers+1)]->setResidualResolution(mwc_resolution);	

		Sensors[(nLayers+2)] = new SensorHitMap((nLayers+2));				//attention: This is specifically tailored for the 8-layer setup
		Sensors[(nLayers+2)]->setLabZ(mwcs->at(1).z, 0.001);
		Sensors[(nLayers+2)]->setCenterHitPosition(mwcs->at(1).x, mwcs->at(1).y ,mwc_resolution , mwc_resolution);
		Sensors[(nLayers+2)]->setParticleEnergy(energy);
		Sensors[(nLayers+2)]->setResidualResolution(mwc_resolution);
	}

	//step 5: fill particle tracks
	Track = new ParticleTrack();
	for (int i=1; i<=nPlanes; i++) {
		Track->addFitPoint(Sensors[i]);
	}
	Track->weightFitPoints(fitPointWeightingMethod);
	Track->fitTrack(fittingMethod);
	
	//step 6: calculate the deviations between each fit missing one layer and exactly that layer's true central position
	if (fittingMethod==GBLTRACK) {
		Track->gblTrackToMilleBinary(milleBinary);
	}else{
		
		for (int layer=1; layer<=nPlanes; layer++) {
			double layer_labZ = Sensors[layer]->getLabZ();
			double intrinsic_z = Sensors[layer]->getIntrinsicHitZPosition();	
			
			std::pair<double, double> position_predicted = Track->calculatePositionXY(layer_labZ+intrinsic_z, layer);
			double x_predicted = position_predicted.first;
			double y_predicted = position_predicted.second;


			std::pair<double, double> position_true = Sensors[layer]->getHitPosition();	
			double x_true = position_true.first;
			double y_true = position_true.second;
			
			Sensors[layer]->getHitPositionError();	

			//step7: calculate the necessary derivatives for Mille and fill them into the binary file			
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
			derGl[(layer-1)*NGLperLayer+2] = y_true;		

			rMeas = x_true - x_predicted;
			sigma = Sensors[layer]->getResidualResolution();
			mille->mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);


			//the y-coordinate
			derLc[0] = 0.;
			derLc[1] = 0.;
			derLc[2] = layer_labZ;
			derLc[3] = 1.;

			derGl[(layer-1)*NGLperLayer+0] = 0.;		
			derGl[(layer-1)*NGLperLayer+1] = 1.;			
			derGl[(layer-1)*NGLperLayer+2] = -x_true;		

			rMeas = y_true - y_predicted;
			sigma = Sensors[layer]->getResidualResolution();
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