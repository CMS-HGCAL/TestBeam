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

double DWC_x0s[4] = {0.04, 0.25, 0.001, 0.0};


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
		
		int nLayers;
		bool useMWCReference;
		bool MWCQualityCut;

		double wc_resolution;
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
		double res1_x, res1_y, res2_x, res2_y, res3_x, res3_y, res4_x, res4_y;
		int NDF;
		double chi2;
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

	nLayers = iConfig.getParameter<int>("nLayers");
	useMWCReference = iConfig.getParameter<bool>("useMWCReference");
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
	if (rd->event>-1)
		return;

	eventCounter++;

	//get the multi wire chambers
	edm::Handle<WireChambers> dwcs;

	event.getByToken(MWCToken, dwcs);
	
	for (int n_layer=0; n_layer<nLayers; n_layer++)
		if (!dwcs->at(n_layer).goodMeasurement) return;
	

	for (size_t n_layer=0; n_layer<4; n_layer++) {
		//Step 4: add dWCs to the setup if useMWCReference option is set true
		Sensors[n_layer] = new SensorHitMap(n_layer);				//attention: This is specifically tailored for the 8-layer setup

		Sensors[n_layer]->setLabZ(dwcs->at(n_layer).z, DWC_x0s[n_layer]);
		Sensors[n_layer]->setCenterHitPosition(dwcs->at(n_layer).x, dwcs->at(n_layer).y , wc_resolution , wc_resolution);
		Sensors[n_layer]->setParticleEnergy(rd->energy);
		Sensors[n_layer]->setResidualResolution(wc_resolution);	
	}

	//first layer of HGCal as reference
	Track = new ParticleTrack();
	Sensors[nLayers] = new SensorHitMap(nLayers);
	Sensors[nLayers]->setLabZ(0., 0.003);
	Sensors[nLayers]->setParticleEnergy(rd->energy);
	Track->addReferenceSensor(Sensors[nLayers]);

	//step 5: fill particle tracks
	for (int n_layer=0; n_layer<nLayers; n_layer++) {
		Track->addFitPoint(Sensors[n_layer]);
	}	
	
	Track->fitTrack(fittingMethod);
	NDF = Track->getNDF();
	chi2 = Track->getChi2();

	std::cout<<"Track chi2: "<<Track->getChi2()<<" / "<<Track->getNDF()<<" = "<<Track->getChi2()/Track->getNDF()<<"   uncertainty at z=0: "<<Track->calculateReferenceErrorXY().first<<" , "<<Track->calculateReferenceErrorXY().second<<std::endl;
	bool accept = true;
	if (MWCQualityCut&&Track->getChi2()/Track->getNDF() > 10) accept=false;		//perform a cut on the quality of the track
	
	if (accept) {
		//step 6: calculate the deviations between each fit missing one layer and exactly that layer's true central position
		for (int n_layer=0; n_layer<nLayers; n_layer++) {
			double layer_labZ = Sensors[n_layer]->getLabZ();
			double intrinsic_z = Sensors[n_layer]->getIntrinsicHitZPosition();	
			
			std::pair<double, double> position_predicted = Track->calculatePositionXY(layer_labZ+intrinsic_z, n_layer);
			double x_predicted = position_predicted.first;
			double y_predicted = position_predicted.second;


			std::pair<double, double> position_true = Sensors[n_layer]->getHitPosition();	
			double x_true = position_true.first;
			double y_true = position_true.second;
			
			if (makeTree) {
				double res_x=x_true - x_predicted;
				double res_y=y_true - y_predicted;
				switch(n_layer) {
					case 0:
						res1_x=res_x;
						res1_y=res_y;
						break;
					case 1:
						res2_x=res_x;
						res2_y=res_y;
						break;
					case 2:
						res3_x=res_x;
						res3_y=res_y;
						break;
					case 3:
						res4_x=res_x;
						res4_y=res_y;
						break;
				}
			}
		}
		tree->Fill();


		if (fittingMethod==GBLTRACK) {
			Track->gblTrackToMilleBinary(milleBinary);
		} else{	
			for (int n_layer=0; n_layer<nLayers; n_layer++) {
				double layer_labZ = Sensors[n_layer]->getLabZ();
				double intrinsic_z = Sensors[n_layer]->getIntrinsicHitZPosition();	
				
				std::pair<double, double> position_predicted = Track->calculatePositionXY(layer_labZ+intrinsic_z, n_layer);
				double x_predicted = position_predicted.first;
				double y_predicted = position_predicted.second;


				std::pair<double, double> position_true = Sensors[n_layer]->getHitPosition();	
				double x_true = position_true.first;
				double y_true = position_true.second;
				

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
				
				//std::cout<<"x_predicted: "<<x_predicted<<"   y_predicted: "<<y_predicted<<"     x_true: "<<x_true<<"   y_true: "<<y_true<<std::endl;

				derGl[n_layer*NGLperLayer+0] = 1.;		
				derGl[n_layer*NGLperLayer+1] = 0.;		
				//derGl[n_layer*NGLperLayer+2] = x_true;		
				//derGl[n_layer*NGLperLayer+3] = 0.;		

				rMeas = x_true - x_predicted;
				//std::cout<<" ---> delta x = "<<rMeas<<"  ";
				sigma = Sensors[n_layer]->getResidualResolution();
				mille->mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);


				//the y-coordinate
				derLc[0] = 0.;
				derLc[1] = 0.;
				derLc[2] = layer_labZ;
				derLc[3] = 1.;

				derGl[n_layer*NGLperLayer+0] = 0.;		
				derGl[n_layer*NGLperLayer+1] = 1.;			
				//derGl[n_layer*NGLperLayer+2] = 0.;			
				//derGl[n_layer*NGLperLayer+3] = y_true;			

				rMeas = y_true - y_predicted;
				//std::cout<<" ---> delta y = "<<rMeas<<std::endl;
				sigma = Sensors[n_layer]->getResidualResolution();
				mille->mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);
			
			}

			mille->end();
		}
	}


	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	delete Track;
	//std::cout<<std::endl;
	
}// analyze ends here


void MillepedeBinaryWriter::beginJob() {	
	wc_resolution=0.5;	//dummy value		//mm

	eventCounter = 0;

	if (fittingMethod==GBLTRACK) {
		milleBinary = new gbl::MilleBinary((binaryFileName).c_str());
		mille = NULL;
	} else{
		milleBinary = NULL;
		mille = new Mille((binaryFileName).c_str());
		std::cout<<"Writing to the file: "<<binaryFileName<<std::endl;
	}
  
	NLC = 4;
	NGLperLayer = 2;	//two translations, no scales, but no rotation
	NGL = nLayers*NGLperLayer;
	rMeas = 0.;
	sigma = 0.;
	derLc = new float[NLC];
	derGl = new float[NGL];
	label = new int[NGL];

	for (int l=0; l<nLayers; l++) {
		label[l*NGLperLayer + 0] = (l+1)*100 + 11;
		label[l*NGLperLayer + 1] = (l+1)*100 + 12;
		//label[l*NGLperLayer + 2] = (l+1)*100 + 21;
		//label[l*NGLperLayer + 3] = (l+1)*100 + 22;
	}


	if (makeTree) {
		tree = fs->make<TTree>("dwc_residuals", "dwc_residuals");
		tree->Branch("res1_x", &res1_x);
		tree->Branch("res1_y", &res1_y);
		tree->Branch("res2_x", &res2_x);
		tree->Branch("res2_y", &res2_y);
		tree->Branch("res3_x", &res3_x);
		tree->Branch("res3_y", &res3_y);
		tree->Branch("res4_x", &res4_x);
		tree->Branch("res4_y", &res4_y);
		tree->Branch("chi2", &chi2);
		tree->Branch("NDF", &NDF);

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