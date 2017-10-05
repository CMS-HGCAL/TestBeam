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
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/Reco/interface/Mille.h"

#include "HGCal/Reco/interface/PositionResolutionHelpers.h"
#include "HGCal/Reco/interface/Tracks.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "Alignment/ReferenceTrajectories/interface/MilleBinary.h"
#include "TTree.h"
#include "TFile.h"

/*
	//DWC E: index 0,
	//DWC D: index 1,
	//DWC A: index 2,
	//DWC ext.: index 3
*/

double DWC_x0s[4] = {0.004, 0.25, 0.008, 0.0};
enum COORDINATE {
	X = 1,
	Y
};


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

		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<WireChambers> MWCToken;
		edm::Service<TFileService> fs;

		int eventCounter;

		std::string methodString;
		TrackFittingMethod fittingMethod;		
		
		std::vector<int> Layers;
		std::string coordinateString;
		COORDINATE coordinate;
		
		bool MWCQualityCut;
		double energy;

		//helper variables that are set within the event loop, i.e. are defined per event
		std::map<int, SensorHitMap*> Sensors;
		ParticleTrack* Track;
		
		Mille* mille;
		gbl::MilleBinary* milleBinary;
		std::string binaryFileName;

  		int NLC, NGLperLayer, NGL;
		float rMeas, sigma;
		float *derLc, *derGl;
		int *label;

		bool makeTree;
		double res1, res2, res3, res4;
		
		TTree* tree;
};

MillepedeBinaryWriter::MillepedeBinaryWriter(const edm::ParameterSet& iConfig) {
	usesResource("TFileService");

	// initialization	
	MWCToken= consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));

	//read the track fitting method
	methodString = iConfig.getParameter<std::string>("fittingMethod");
	if (methodString == "lineAnalytical")
		fittingMethod = LINEFITANALYTICAL;
	else 
		fittingMethod = DEFAULTFITTING;

	Layers = iConfig.getParameter<std::vector<int> >("Layers");
	
	coordinateString = iConfig.getParameter<std::string>("Coordinate");
	std::cout<<"CoordinateString: "<<coordinateString<<std::endl;
	if (coordinateString == "x") coordinate = X;
	else if (coordinateString == "y") coordinate = Y;
	else coordinate = X;

	MWCQualityCut = iConfig.getParameter<bool>("MWCQualityCut");

	binaryFileName = iConfig.getParameter<std::string>("binaryFile");

	makeTree = iConfig.getUntrackedParameter<bool>("makeTree", true);
	tree = NULL;


}//constructor ends here

MillepedeBinaryWriter::~MillepedeBinaryWriter() {
	return;
}

// ------------ method called for each event  ------------
void MillepedeBinaryWriter::analyze(const edm::Event& event, const edm::EventSetup& setup) {
	edm::Handle<RunData> rd;
 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	
	eventCounter++;

	//get the multi wire chambers
	edm::Handle<WireChambers> dwcs;
	event.getByToken(MWCToken, dwcs);
	

	for (int k=0; k<NGL; k++){
		derGl[k] = 0.;
	}
	res1 = res2 = res3 = res4 = -999;


	//Step 4: add dWCs to the setup if useMWCReference option is set true
	int N_points = 0;
	Track = new ParticleTrack();
	Sensors[100] = new SensorHitMap(100);
	Sensors[100]->setLabZ(0., 0.003);
	Sensors[100]->setParticleEnergy(rd->energy);
	Track->addReferenceSensor(Sensors[100]);
	
	//first the 
	for (size_t l=0; l<Layers.size(); l++) {
		int n_layer = Layers[l];

		if (!dwcs->at(n_layer).goodMeasurement) return;
		
		Sensors[n_layer] = new SensorHitMap(n_layer);				
		Sensors[n_layer]->setLabZ(dwcs->at(n_layer).z, DWC_x0s[n_layer]);
		Sensors[n_layer]->setCenterHitPosition(dwcs->at(n_layer).x, dwcs->at(n_layer).y , dwcs->at(n_layer).res_x , dwcs->at(n_layer).res_y);
		Sensors[n_layer]->setParticleEnergy(rd->energy);
		//Sensors[n_layer]->setResidualResolution(dwcs->at(n_layer).res_x);	
		Sensors[n_layer]->setResidualResolution(3.);	

		Track->addFitPoint(Sensors[n_layer]);
		N_points++;
	}
	
	//quality cut to remove unphysical events:
	double _c0 = coordinate==X ? Sensors[0]->getHitPosition().first : Sensors[0]->getHitPosition().second;
	double _c1 = coordinate==X ? Sensors[1]->getHitPosition().first : Sensors[1]->getHitPosition().second;
	double _c2 = coordinate==X ? Sensors[2]->getHitPosition().first : Sensors[2]->getHitPosition().second;
	double _c3 = coordinate==X ? Sensors[3]->getHitPosition().first : Sensors[3]->getHitPosition().second;
	if (fabs(_c0-_c1)>10.) return;
	if (fabs(_c0-_c2)>50.) return;
	if (fabs(_c0-_c3)>50.) return;
	if (fabs(_c1-_c2)>50.) return;
	if (fabs(_c1-_c3)>50.) return;
	if (fabs(_c2-_c3)>10.) return;


	if (N_points == (int)Layers.size()) {
		Track->fitTrack(fittingMethod);
		
		if (MWCQualityCut&&fabs(Track->getChi2(coordinate)/Track->getNDF(coordinate)) > 30.) return;
	
		//step 6: calculate the deviations between each fit missing one layer and exactly that layer's true central position
		for (size_t l=0; l<Layers.size(); l++) {
			int n_layer = Layers[l];
			
			double layer_labZ = Sensors[n_layer]->getLabZ();
			double intrinsic_z = Sensors[n_layer]->getIntrinsicHitZPosition();	
			
			double _predicted = coordinate==X ? Track->calculatePositionXY(layer_labZ+intrinsic_z, n_layer).first : Track->calculatePositionXY(layer_labZ+intrinsic_z, n_layer).second;
			double _true = coordinate==X ? Sensors[n_layer]->getHitPosition().first : Sensors[n_layer]->getHitPosition().second;
			
			//the x-coordinate								
			derLc[0] = layer_labZ;
			derLc[1] = 1.;
			
			//std::cout<<"x_predicted: "<<x_predicted<<"   y_predicted: "<<y_predicted<<"     x_true: "<<x_true<<"   y_true: "<<y_true<<std::endl;
			for (int k=0; k<NGL; k++){
				derGl[k] = 0.;
			}

			derGl[n_layer*NGLperLayer+0] = 1.;		
				
			rMeas = _true - _predicted;
			
			//std::cout<<" ---> delta x = "<<rMeas<<"  ";
			sigma = Sensors[n_layer]->getResidualResolution();
			mille->mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);
			//std::cout<<"Adding rMeas (in x) = "<<rMeas<<std::endl;

			if (makeTree) {
				double res=_true - _predicted;
				switch(n_layer) {
					case 0:
						res1=res;
						if (fabs(res1)>30.) {
							std::cout<<"Layer 1: "<<Sensors[0]->getHitPosition().first<<"   "<<Sensors[0]->getHitPosition().second<<std::endl;
							std::cout<<"Layer 2: "<<Sensors[1]->getHitPosition().first<<"   "<<Sensors[1]->getHitPosition().second<<std::endl;
							std::cout<<"Layer 3: "<<Sensors[2]->getHitPosition().first<<"   "<<Sensors[2]->getHitPosition().second<<std::endl;
							std::cout<<"Layer 4: "<<Sensors[3]->getHitPosition().first<<"   "<<Sensors[3]->getHitPosition().second<<std::endl;
						}
						break;
					case 1:
						res2=res;
						break;
					case 2:
						res3=res;
						break;
					case 3:
						res4=res;
						break;
				}
			}
		}
		mille->end();
		tree->Fill();	
	}
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	delete Track;
		

	
}// analyze ends here


void MillepedeBinaryWriter::beginJob() {	

	eventCounter = 0;

	if (fittingMethod==GBLTRACK) {
		milleBinary = new gbl::MilleBinary((binaryFileName).c_str());
		mille = NULL;
	} else{
		milleBinary = NULL;
		mille = new Mille((binaryFileName).c_str());
		std::cout<<"Writing to the file: "<<binaryFileName<<std::endl;
	}
  
	NLC = 2;
	NGLperLayer = 1;	//two translations, no scales, but no rotation
	NGL = 4*NGLperLayer;		//hard coded number for four DWCs
	std::cout<<"Number of global parameters: "<<NGL<<std::endl;
	rMeas = 0.;
	sigma = 0.;
	derLc = new float[NLC];
	derGl = new float[NGL];
	label = new int[NGL];

	for (size_t l=0; l<Layers.size(); l++) {
		int layer = Layers[l];
		label[layer*NGLperLayer + 0] = (layer)*100 + (coordinate==X ? 11 : 12);
	}

	if (makeTree) {
		tree = fs->make<TTree>("dwc_residuals", "dwc_residuals");
		tree->Branch((std::string("res1_")+std::string(coordinate==X?"x":"y")).c_str(), &res1);
		tree->Branch((std::string("res2_")+std::string(coordinate==X?"x":"y")).c_str(), &res2);
		tree->Branch((std::string("res3_")+std::string(coordinate==X?"x":"y")).c_str(), &res3);
		tree->Branch((std::string("res4_")+std::string(coordinate==X?"x":"y")).c_str(), &res4);
	}
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
	
}

void MillepedeBinaryWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MillepedeBinaryWriter);