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
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Tracks.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "TFile.h"
#include "TTree.h"
  
//configuration1 based on Shilpis plot from 26.07.17: 
//https://indico.cern.ch/event/656159/contributions/2674177/attachments/1499439/2334623/update_simulation_geom.pdf
// 2 layers in EE (first one in the graphics was removed) and 4 in FH, indications are in mm
double z1_config1 = 0. + 156. + 182.7 + 4.9 * 4 + 6.0 * 4 + 0.3 * 7 + 1.2;
double z2_config1 = z1_config1 + 11.0 + 4.9 * 8 + 0.3 * 16 + 6.0 * 7 + 13.2;
double z3_config1 = z2_config1 + 1.2 + 6.0 * 3 + 4 * 0.3 + 2 * 4.9 + 54. + 38. + 10. + 78. + 40. + 6. + 1.2;	
double z4_config1 = z3_config1 + 46.9 + 3 * 40 + 6 + 1.2;
double z5_config1 = z4_config1 + 26.9 + 40. + 10. + 63. + 10. + 2 * 40. + 6. + 1.2;
double z6_config1 = z5_config1 + 58.9 + 3 * 40. + 6. + 1.2;

double zdepth1_config1 = 0.;
double zdepth2_config1 = 0.; 
double zdepth3_config1 = 0.; 
double zdepth4_config1 = 0.; 
double zdepth5_config1 = 0.; 
double zdepth6_config1 = 0.; 

//here they must be indicated in cm
double config1Positions[] = {z1_config1/10., z2_config1/10., z3_config1/10., z4_config1/10., z5_config1/10., z6_config1/10.};    
double  config1X0Depths[] = {zdepth1_config1, zdepth2_config1, zdepth3_config1, zdepth4_config1, zdepth5_config1, zdepth6_config1}; 
  


//#define DEBUG

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
	 	
		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<WireChambers> MWCToken;
		
		AlignmentParameters* alignmentParameters; //all entries are set to zero if no valid file is given 

		ConsiderationMethod considerationMethod;
		WeightingMethod weightingMethod;
		TrackFittingMethod fittingMethod;		

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
		int configuration, evId, eventCounter, run, pdgID, layer; 	//eventCounter: counts the events in this analysis run to match information within ove event to each other
		double energy;
		double layerWeight, layerEnergy, layerClusterEnergy, sumFitWeights, sumEnergy, sumClusterEnergy, CM_cells_count, CM_sum;
		double track_chi2_x, track_chi2_y;
		int track_ndf_x, track_ndf_y;
		int NgoodDWCs, DWCReferenceType;
		double x_predicted, x_predicted_err, y_predicted, y_predicted_err, x_true, x_true_err, y_true, y_true_err, deltaX, deltaY;
		double x_predicted_to_closest_cell, y_predicted_to_closest_cell, x_true_to_closest_cell, y_true_to_closest_cell, layerZ_cm, layerZ_X0, deviation;

		//averaged information up to the corresponding layers
		double average_x_predicted, average_y_predicted, average_x_true, average_y_true, average_deltaX, average_deltaY; 

		//dWCs
		int useMWC;

		std::pair<int, double> CM_tmp;	//will write the subtract_CM() return values for each layer
};

Position_Resolution_Analyzer::Position_Resolution_Analyzer(const edm::ParameterSet& iConfig) {	
	
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	MWCToken= consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));

	//read the cell consideration option to calculate the central hit point
	std::string methodString = iConfig.getParameter<std::string>("considerationMethod");
	if (methodString == "all")
		considerationMethod = CONSIDERALL;
	else if (methodString == "closest7")
		considerationMethod = CONSIDERSEVEN;
	else if (methodString == "closest19")
		considerationMethod = CONSIDERNINETEEN;

	//read the weighting method to obtain the central hit point
	methodString = iConfig.getParameter<std::string>("weightingMethod");
	if (methodString == "squaredWeighting")
		weightingMethod = SQUAREDWEIGHTING;	
	else if (methodString == "linearWeighting")
		weightingMethod = LINEARWEIGHTING;
	else if (methodString == "logWeighting_3.5_1.0")
		weightingMethod = LOGWEIGHTING_35_10;
	else if (methodString == "mostIntensive")
		weightingMethod = MOSTINTENSIVE;
	else 
		weightingMethod = DEFAULTWEIGHTING;

	//read the track fitting method
	methodString = iConfig.getParameter<std::string>("fittingMethod");
	if (methodString == "lineAnalytical")
		fittingMethod = LINEFITANALYTICAL;
	else if (methodString == "lineTGraphErrors")
		fittingMethod = LINEFITTGRAPHERRORS;
	else if (methodString == "gblTrack")
		fittingMethod = GBLTRACK;
	else 
		fittingMethod = DEFAULTFITTING;


	//read the layer configuration
	LayersConfig = iConfig.getParameter<int>("layers_config");
	if (LayersConfig == 1) {
		Layer_Z_Positions = std::vector<double>(config1Positions, config1Positions + sizeof(config1Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(config1X0Depths, config1X0Depths + sizeof(config1X0Depths)/sizeof(double));
	} else {
		Layer_Z_Positions = std::vector<double>(config1Positions, config1Positions + sizeof(config1Positions)/sizeof(double));
		Layer_Z_X0s 			= std::vector<double>(config1X0Depths, config1X0Depths + sizeof(config1X0Depths)/sizeof(double));
	}

	eventCounter = 0;

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
	outTree->Branch("pdgID", &pdgID, "pdgID/I");
	outTree->Branch("layer", &layer, "layer/I");
	outTree->Branch("energy", &energy, "energy/D");	//electron energy in GeV
	
	outTree->Branch("layerEnergy", &layerEnergy, "layerEnergy/D");
	outTree->Branch("sumEnergy", &sumEnergy, "sumEnergy/D");
	
	outTree->Branch("x_true", &x_true, "x_true/D");
	outTree->Branch("x_true_to_closest_cell", &x_true_to_closest_cell, "x_true_to_closest_cell/D");
	outTree->Branch("x_true_err", &x_true_err, "x_true_err/D");
	outTree->Branch("y_true", &y_true, "y_true/D");
	outTree->Branch("y_true_to_closest_cell", &y_true_to_closest_cell, "y_true_to_closest_cell/D");
	outTree->Branch("y_true_err", &y_true_err, "y_true_err/D");
	
	outTree->Branch("NgoodDWCs", &NgoodDWCs, "NgoodDWCs/I");
	outTree->Branch("DWCReferenceType", &DWCReferenceType, "DWCReferenceType/I");
	outTree->Branch("track_chi2_x", &track_chi2_x, "track_chi2_x/D");
	outTree->Branch("track_chi2_y", &track_chi2_y, "track_chi2_y/D");
	outTree->Branch("track_ndf_x", &track_ndf_x, "track_ndf_x/I");
	outTree->Branch("track_ndf_y", &track_ndf_y, "track_ndf_y/I");

	outTree->Branch("x_predicted", &x_predicted, "x_predicted/D");
	outTree->Branch("x_predicted_to_closest_cell", &x_predicted_to_closest_cell, "x_predicted_to_closest_cell/D");
	outTree->Branch("x_predicted_err", &x_predicted_err, "x_predicted_err/D");
	outTree->Branch("y_predicted", &y_predicted, "y_predicted/D");
	outTree->Branch("y_predicted_to_closest_cell", &y_predicted_to_closest_cell, "y_predicted_to_closest_cell/D");
	outTree->Branch("y_predicted_err", &y_predicted_err, "y_predicted_err/D");
	
	outTree->Branch("deltaX", &deltaX, "deltaX/D");
	outTree->Branch("deltaY", &deltaY, "deltaY/D");
	outTree->Branch("deviation", &deviation, "deviation/D");

	outTree->Branch("average_x_predicted", &average_x_predicted, "average_x_predicted/D");
	outTree->Branch("average_y_predicted", &average_y_predicted, "average_y_predicted/D");
	outTree->Branch("average_x_true", &average_x_true, "average_x_true/D");
	outTree->Branch("average_y_true", &average_y_true, "average_y_true/D");
	outTree->Branch("average_deltaX", &average_deltaX, "average_deltaX/D");
	outTree->Branch("average_deltaY", &average_deltaY, "average_deltaY/D");

	alignmentParameters = new AlignmentParameters(iConfig.getParameter<std::vector<std::string> >("alignmentParameterFiles")); 

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
	pdgID = std::atoi( (rd->runType).c_str() );
	eventCounter = rd->event;
	energy = rd->energy;
	
	if (rd->hasDanger) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because somthing went wrong"<<std::endl;
		return;
	}

	if (run == -1) {
		std::cout<<"Run is not in configuration file - is ignored."<<std::endl;
		return;
	}

	edm::Handle<WireChambers> dwcs;
	event.getByToken(MWCToken, dwcs);

	//initialize new fit counters in case this is a new run:
	if (successfulFitCounter.find(run) == successfulFitCounter.end()) 
		successfulFitCounter[run] = failedFitCounter[run] = 0;

	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);


	//step 1: Reduce the information to energy deposits/hits in x,y per sensor/layer 
	//fill the rechits:
	for(auto Rechit : *Rechits) {	
		layer = (Rechit.id()).layer();
  
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap(layer);
			Sensors[layer]->setLabZ(Layer_Z_Positions[layer-1], Layer_Z_X0s[layer-1]);	//first argument: real positon as measured (not aligned) in cm, second argument: position in radiation lengths

			Sensors[layer]->setAlignmentParameters(alignmentParameters->getValue(run, 100*layer + 21), 0.0, 0.0,
				alignmentParameters->getValue(run, 100*layer + 11), alignmentParameters->getValue(run, 100*layer + 12), 0.0);	
			Sensors[layer]->setSensorSize(SensorSize);

			double X0sum = 0;
			for (int _x = 0; _x<(int)layer; _x++) X0sum += Layer_Z_X0s[_x];
			Sensors[layer]->setParticleEnergy(energy - gblhelpers::computeEnergyLoss(X0sum, energy));
		}

		Sensors[layer]->addHit(Rechit, 1.0);		
	}

	#ifdef DEBUG
		std::cout<<"run: "<<rd->run<<"  energy: "<<rd->energy<<"  type:" << rd->runType<<"   eventCounter: "<<rd->event<<std::endl;
	#endif
	
	//Possible event selection: sum of energies of all cells(=hits) from RecHits Collection and Clusters
	sumEnergy = 0.;
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {	
		sumEnergy += it->second->getTotalEnergy();
	}

	

	//step 2: calculate impact point with technique indicated as the argument
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(considerationMethod, weightingMethod);
		
		#ifdef DEBUG
			std::pair<double, double> position_true = it->second->getLabHitPosition();
			std::cout<<"layer "<<it->first<<"  x: "<<position_true.first<<"    y: "<<position_true.second<<std::endl;
		#endif
	}

	#ifdef DEBUG
		std::cout<<dwcs->at(i).x<<"  "<<dwcs->at(0).y<<"   "<<dwcs->at(0).z<<"  "<<dwcs->at(0).goodMeasurement<<std::endl;
		std::cout<<dwcs->at(1).x<<"  "<<dwcs->at(1).y<<"   "<<dwcs->at(1).z<<"  "<<dwcs->at(1).goodMeasurement<<std::endl;
		std::cout<<dwcs->at(2).x<<"  "<<dwcs->at(2).y<<"   "<<dwcs->at(2).z<<"  "<<dwcs->at(2).goodMeasurement<<std::endl;
		std::cout<<dwcs->at(3).x<<"  "<<dwcs->at(3).y<<"   "<<dwcs->at(3).z<<"  "<<dwcs->at(3).goodMeasurement<<std::endl;
	#endif

	//Step 3: add dWCs to the setup if useMWCReference option is set true
	if (useMWCReference) {
		useMWC = 1;
		NgoodDWCs = 0;
		DWCReferenceType = 0;
		

		for (size_t i=0; i<dwcs->size(); i++) {
			if (!dwcs->at(i).goodMeasurement) continue;
			
			NgoodDWCs++;
			DWCReferenceType += pow(2, i);
			
			Sensors[(nLayers+NgoodDWCs)] = new SensorHitMap((nLayers+NgoodDWCs));				//attention: This is specifically tailored for the 8-layer setup
			Sensors[(nLayers+NgoodDWCs)]->setLabZ(dwcs->at(i).z, 0.001);
			Sensors[(nLayers+NgoodDWCs)]->setCenterHitPosition(dwcs->at(i).x, dwcs->at(i).y, dwcs->at(i).res_x , dwcs->at(i).res_y);
			Sensors[(nLayers+NgoodDWCs)]->setParticleEnergy(energy);
			Sensors[(nLayers+NgoodDWCs)]->setAlignmentParameters(alignmentParameters->getValue(energy, 100*(nLayers+NgoodDWCs) + 21), 0.0, 0.0,
					alignmentParameters->getValue(energy, 100*(nLayers+NgoodDWCs) + 11), alignmentParameters->getValue(energy, 100*(nLayers+NgoodDWCs) + 12), 0.0);	
			Sensors[(nLayers+NgoodDWCs)]->setResidualResolution(dwcs->at(i).res_x);	
		}

		#ifdef DEBUG
			std::cout<<"DWC reference type: "<<DWCReferenceType<<std::endl;
		#endif
		if (NgoodDWCs<2) return;

	} else {
		useMWC = 0;
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
				if(useMWCReference && j>nLayers) {
					//std::cout<<"Adding sensor "<<j<<" to track i"<<std::endl;
					Tracks[i]->addFitPoint(Sensors[j]);
				}
				else if(!useMWCReference && j<=nLayers) Tracks[i]->addFitPoint(Sensors[j]);
			} else {
				if(useMWCReference && j<=nLayers) 
					Tracks[i]->addFitPoint(Sensors[j]);					
			}
		}

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
		x_predicted = -position_predicted.first;
		sum_x_predicted  += x_predicted*layerEnergy;
		y_predicted = position_predicted.second;
		sum_y_predicted  += y_predicted*layerEnergy;

		track_chi2_x = Tracks[layer]->getChi2(1);
		track_chi2_y = Tracks[layer]->getChi2(2);
		track_ndf_x = Tracks[layer]->getNDF(1);
		track_ndf_y = Tracks[layer]->getNDF(2);

		std::pair<double, double> position_predicted_to_closest_cell = Sensors[layer]->getCenterOfClosestCell(position_predicted);
		x_predicted_to_closest_cell = -position_predicted_to_closest_cell.first;	//2017: DWCs are rotated by 180 degrees in x
		y_predicted_to_closest_cell = position_predicted_to_closest_cell.second;
		std::pair<double, double> position_error_predicted = Tracks[layer]->calculateReferenceErrorXY();
		x_predicted_err = position_error_predicted.first;		//2017: DWCs are rotated by 180 degrees in x
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

		average_x_predicted = layer <= 6 ? sum_x_predicted / sum_energy : -999;
		average_y_predicted = layer <= 6 ? sum_y_predicted / sum_energy : -999;
		average_x_true = layer <= 6 ? sum_x_true / sum_energy : -999;
		average_y_true = layer <= 6 ? sum_y_true / sum_energy : -999;
		average_deltaX = layer <= 6 ? average_x_true - average_x_predicted : -999;
		average_deltaY = layer <= 6 ? average_y_true - average_y_predicted : -999;
		#ifdef DEBUG
			if (layer==1) {
				std::cout<<"   layer: "<<layer<<"   x:  "<<x_predicted<<" - "<<x_true<<"     "<<y_predicted<<" - "<<y_true<<std::endl;
			}
		#endif
		
		outTree->Fill();
		//fill the tree
	}
	/*
	*/
	
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	//for(std::map<int, ParticleTrack*>::iterator it=Tracks.begin(); it!=Tracks.end(); it++) {
	//	delete (*it).second;
	//}; Tracks.clear();
	
}// analyze ends here

void Position_Resolution_Analyzer::beginJob() {	
}

void Position_Resolution_Analyzer::endJob() {
	

	delete alignmentParameters;
}

void Position_Resolution_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Position_Resolution_Analyzer);