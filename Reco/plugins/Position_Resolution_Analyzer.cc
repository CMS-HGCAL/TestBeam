/* Write  
 * a description here!
 */

/**
	@Author: Thorben Quast <tquast>
		16 Nov 2016
		thorben.quast@cern.ch / thorben.quast@rwth-aachen.de
*/



// system include files
//#include <memory>
#include <iostream>
#include <string>
#include <vector>
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
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph2D.h"
//#include "TFileDirectory.h"



/**********/
//all of this should go into a config file

double Layer_Z_Positions[17]  = {0.0, 1.2, 2., 3.5, 4.3, 5.8, 6.3, 8.7, 9.5, 11.4, 12.2, 13.8, 14.6, 16.6, 17.4, 20., 20.8};
//by this convention the first entry is not considered as the layers numbering starts at 1
int nLayers = 8;
int SensorSize = 128;
/**********/
                    
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
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_;
		edm::EDGetTokenT<RunData> RunDataToken;	
		
		WeightingMethod weightingMethod;
		TrackFittingMethod fittingMethod;		
		bool make2DGraphs;

		int successfulFitCounter, failedFitCounter;
		double min_deviation, max_deviation; 	//minimum and maximum value of the deviations for the 2D histogram, those are defined globally

		//stuff to be written to objects at the end, first is the energy, second the layer thickness, third the run and last the layer int
		std::map<double, std::map<double, std::map<int, std::map<int, std::vector<double> > > > >deviations;

		//helper variables that are set within the event loop, i.e. are defined per event
		int evId, run, layer;
		double energy, layerThickness;
		std::map<int, SensorHitMap*> Sensors;
		std::map<int, ParticleTrack*> Tracks;
		std::vector<double> x_predicted_v,y_predicted_v, x_true_v, y_true_v, layerZ_v;
		double x_predicted, y_predicted, x_true, y_true, layerZ, deviation;
};

Position_Resolution_Analyzer::Position_Resolution_Analyzer(const edm::ParameterSet& iConfig) {
	gStyle->SetOptStat();
	// initialization
	usesResource("TFileService");
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));

	//read the weighting method to obtain the central hit point
	std::string methodString = iConfig.getParameter<std::string>("weightingMethod");
	if (methodString == "squaredWeighting")
		weightingMethod = SQUAREDWEIGHTING;	
	else if (methodString == "linearWeighting")
		weightingMethod = LINEARWEIGHTING;
	else 
		weightingMethod = DEFAULTWEIGHTING;

	//read the track fitting method
	methodString = iConfig.getParameter<std::string>("fittingMethod");
	if (methodString == "lineTGraphErrors")
		fittingMethod = LINEFITTGRAPHERRORS;
	else 
		fittingMethod = DEFAULTFITTING;

	//making 2DGraphs per event?
	make2DGraphs = iConfig.getParameter<bool>("make2DGraphs");

	//initiate some counters that are printed at the end
	successfulFitCounter = failedFitCounter = 0;

	//initiate the minimum and maximum values for the deviation
	min_deviation = pow(10., 12);
	max_deviation = -1.;


}//constructor ends here

Position_Resolution_Analyzer::~Position_Resolution_Analyzer() {
	return;
}

// ------------ method called for each event  ------------
void Position_Resolution_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {

 	//get the relevant event information
	evId = event.id().event();
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);
	run = rd->run;
	energy = rd->energy;
	layerThickness = rd->layerThickness;
	if (run == -1) {
		std::cout<<"Run is not in configuration file - is ignored."<<std::endl;
		return;
	}

	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_, Rechits);
	
	//step 1: Reduce the information to energy deposits/hits in x,y per sensor/layer 
	for(auto Rechit : *Rechits) {	
		layer = (Rechit.id()).layer();
		if ( Sensors.find(layer) == Sensors.end() ) {
			Sensors[layer] = new SensorHitMap();
			Sensors[layer]->setZ(Layer_Z_Positions[layer]);
			Sensors[layer]->setSensorSize(SensorSize);
		}
		Sensors[layer]->addHit(Rechit);
	}

	//step 2: calculate impact point with technique indicated as the argument
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		//subtract pedestals first
		it->second->subtractPedestals();

		//now calculate the center positions for each layer
		it->second->calculateCenterPosition(weightingMethod);
		
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
	for (int i=1; i<=nLayers; i++) {
		layerZ = Sensors[i]->getZ();
		x_predicted = Tracks[i]->calculatePositionXY(layerZ).first;
		y_predicted = Tracks[i]->calculatePositionXY(layerZ).second;
		if (x_predicted==0 && y_predicted==0)	{
			//default fitting has been applied, i.e. the regular fit has failed or the selected method is not implemented
			failedFitCounter++;
			continue; 	//ignore those cases but count them
		}
		successfulFitCounter++; 
		
		x_true = Sensors[i]->getCenterPosition().first;
		y_true = Sensors[i]->getCenterPosition().second;

		deviation  = pow(x_predicted - x_true, 2); 
		deviation += pow(y_predicted - y_true, 2);
		deviation  = sqrt(deviation);
		deviations[energy][layerThickness][run][i].push_back(deviation);

		//update the deviations on the fly
		min_deviation = min_deviation > deviation ? deviation: min_deviation;
		max_deviation = max_deviation < deviation ? deviation: max_deviation;
	
		//store for the two 2D graphs that are written per event
		x_predicted_v.push_back(x_predicted);
		y_predicted_v.push_back(y_predicted);
		x_true_v.push_back(x_true);
		y_true_v.push_back(y_true);
		layerZ_v.push_back(layerZ);
	}

	if (make2DGraphs) {
		std::string graphIdentifier = "run_" + std::to_string(run) + "event_" + std::to_string(evId);
		fs->make<TGraph2D>(("predicted_points_" + graphIdentifier).c_str(), "", layerZ_v.size(), &(x_predicted_v[0]), &(y_predicted_v[0]), &(layerZ_v[0]));
		fs->make<TGraph2D>(("true_points_" + graphIdentifier).c_str(), "", layerZ_v.size(), &(x_true_v[0]), &(y_true_v[0]), &(layerZ_v[0]));
	}

	//Todo: 1. cleanup to avoid memory leak for many events, try to install cmssw on local machine for that!!!
	//2. look at Andre's weighting/fitting

	//clear:
	//loop over entries and delete manually by hand before clearing the maps
	
	
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	
	
	for(std::map<int, ParticleTrack*>::iterator it=Tracks.begin(); it!=Tracks.end(); it++) {
		delete (*it).second;
	}; Tracks.clear();
	
	x_predicted_v.clear(); y_predicted_v.clear(); x_true_v.clear(); y_true_v.clear(); layerZ_v.clear();

}// analyze ends here

void Position_Resolution_Analyzer::beginJob() {	
}

void Position_Resolution_Analyzer::endJob() {
	std::cout<<"*************************************************"<<std::endl;
	std::cout<<"END OF FITTING:"<<std::endl<<std::endl;
	std::cout<<"Successful fits: "<<successfulFitCounter<<std::endl;
	std::cout<<"Failed fits: "<<failedFitCounter<<std::endl;
	
	std::cout<<"Making deviation TH2D(s)... "<<std::endl;
	
	//std::map<double, std::map<double, std::map<int, std::vector<double> > > >deviations;
	std::map<double, std::map<double, std::map<int, std::map<int, std::vector<double> > > > >::iterator it1;
	TFileDirectory subDir1;
	std::map<double, std::map<int, std::map<int, std::vector<double> > > >::iterator it2;
	TFileDirectory subDir2;
	std::map<int, std::map<int, std::vector<double> > >::iterator it3;

	for (it1=deviations.begin(); it1!=deviations.end(); it1++) {
		subDir1 = this->fs->mkdir(std::to_string((*it1).first).c_str());
		for (it2=(*it1).second.begin(); it2!=(*it1).second.end(); it2++) {
			subDir2 = subDir1.mkdir(std::to_string((*it2).first).c_str());
			for (it3=(*it2).second.begin(); it3!=(*it2).second.end(); it3++) {
				TH2D* deviationHistogram = subDir2.make<TH2D>(("run_"+std::to_string((*it3).first)).c_str(), "", nLayers, 0.5, nLayers+0.5, 1000, min_deviation, max_deviation);
				deviationHistogram->GetXaxis()->SetTitle("n_{Layer}");
				deviationHistogram->GetYaxis()->SetTitle("deviation_{x-y} [cm]");
				for(int i=1; i<=nLayers; i++)
					for(int j=0; j<(int)(*it3).second[i].size(); j++)
						deviationHistogram->Fill(i, (*it3).second[i][j]);
			}
		}
	}
	std::cout<<"*************************************************"<<std::endl;
}

void Position_Resolution_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Position_Resolution_Analyzer);