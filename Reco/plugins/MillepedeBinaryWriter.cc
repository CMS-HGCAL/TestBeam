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
		int acceptedCounter;

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
		double reco1_x, reco1_y, reco2_x, reco2_y, reco3_x, reco3_y, reco4_x, reco4_y;
	  	double time_DWC1, time_DWC2, time_DWC3, time_DWC4;
	  	double x1_m_x2, x1_m_x3, x1_m_x4, x2_m_x3, x2_m_x4, x3_m_x4;
	  	double y1_m_y2, y1_m_y3, y1_m_y4, y2_m_y3, y2_m_y4, y3_m_y4;
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

	eventCounter++;

	//get the multi wire chambers
	edm::Handle<WireChambers> dwcs;

	event.getByToken(MWCToken, dwcs);
	
	for (int n_layer=0; n_layer<nLayers; n_layer++)
		if (!dwcs->at(n_layer).goodMeasurement) return;
	

	
	for (size_t n_layer=0; n_layer<4; n_layer++) {
		//Step 4: add dWCs to the setup if useMWCReference option is set true
		Sensors[n_layer] = new SensorHitMap(n_layer);				//attention: This is specifically tailored for the 8-layer setup

		Sensors[n_layer]->setLabZ(dwcs->at(n_layer).z, 1.);
		Sensors[n_layer]->setCenterHitPosition(dwcs->at(n_layer).x, dwcs->at(n_layer).y ,wc_resolution , wc_resolution);
		Sensors[n_layer]->setParticleEnergy(energy);
		Sensors[n_layer]->setResidualResolution(wc_resolution);	
	}

	//step 5: fill particle tracks
	Track = new ParticleTrack();
	for (int n_layer=0; n_layer<nLayers; n_layer++) {
		Track->addFitPoint(Sensors[n_layer]);
	}	
	
	Track->fitTrack(fittingMethod);

	
	//step 6: calculate the deviations between each fit missing one layer and exactly that layer's true central position
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
			
			if (MWCQualityCut&&(fabs(x_true - x_predicted) > 10. || fabs(y_true - y_predicted) > 10.)) return;

			Sensors[n_layer]->getHitPositionError();	

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

		if (makeTree) {

		 	time_DWC1 = dwcs->at(0).averagedTimeStamp;
		 	time_DWC2 = dwcs->at(1).averagedTimeStamp;
		 	time_DWC3 = dwcs->at(2).averagedTimeStamp;
		 	time_DWC4 = dwcs->at(3).averagedTimeStamp;

		 	x1_m_x2 = dwcs->at(0).x - dwcs->at(1).x;
		 	x1_m_x3 = dwcs->at(0).x - dwcs->at(2).x;
			x1_m_x4 = dwcs->at(0).x - dwcs->at(3).x;
			x2_m_x3 = dwcs->at(1).x - dwcs->at(2).x;
			x2_m_x4 = dwcs->at(1).x - dwcs->at(3).x;
			x3_m_x4 = dwcs->at(2).x - dwcs->at(3).x;
			y1_m_y2 = dwcs->at(0).y - dwcs->at(1).y;
			y1_m_y3 = dwcs->at(0).y - dwcs->at(2).y;
			y1_m_y4 = dwcs->at(0).y - dwcs->at(3).y;
			y2_m_y3 = dwcs->at(1).y - dwcs->at(2).y;
			y2_m_y4 = dwcs->at(1).y - dwcs->at(3).y;
			y3_m_y4 = dwcs->at(2).y - dwcs->at(3).y;
	 
			reco1_x = dwcs->at(0).x;
			reco1_y = dwcs->at(0).y;
			reco2_x = dwcs->at(1).x;
			reco2_y = dwcs->at(1).y;
			reco3_x = dwcs->at(2).x;
			reco3_y = dwcs->at(2).y;
			reco4_x = dwcs->at(3).x;
			reco4_y = dwcs->at(3).y;

			tree->Fill();
		}
		mille->end();
	}
	
	for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
		delete (*it).second;
	};	Sensors.clear();
	delete Track;
	//std::cout<<std::endl;
	
	if ((rd->event!=-1)&&(rd->runType==1))
		acceptedCounter++;
}// analyze ends here


void MillepedeBinaryWriter::beginJob() {	
	wc_resolution=1.0;
	energy=250;

	acceptedCounter = 0;
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
		tree = fs->make<TTree>("dwc_reco", "dwc_reco");
		tree->Branch("res1_x", &res1_x);
		tree->Branch("res1_y", &res1_y);
		tree->Branch("res2_x", &res2_x);
		tree->Branch("res2_y", &res2_y);
		tree->Branch("res3_x", &res3_x);
		tree->Branch("res3_y", &res3_y);
		tree->Branch("res4_x", &res4_x);
		tree->Branch("res4_y", &res4_y);
		tree->Branch("reco1_x", &reco1_x);
		tree->Branch("reco1_y", &reco1_y);
		tree->Branch("reco2_x", &reco2_x);
		tree->Branch("reco2_y", &reco2_y);
		tree->Branch("reco3_x", &reco3_x);
		tree->Branch("reco3_y", &reco3_y);
		tree->Branch("reco4_x", &reco4_x);
		tree->Branch("reco4_y", &reco4_y);
	  	tree->Branch("time_DWC1", &time_DWC1);
	  	tree->Branch("time_DWC2", &time_DWC2);
	  	tree->Branch("time_DWC3", &time_DWC3);
	  	tree->Branch("time_DWC4", &time_DWC4);
	  	tree->Branch("x1_m_x2", &x1_m_x2);
	  	tree->Branch("x1_m_x3", &x1_m_x3);
	  	tree->Branch("x1_m_x4", &x1_m_x4);
	  	tree->Branch("x2_m_x3", &x2_m_x3);
	  	tree->Branch("x2_m_x4", &x2_m_x4);
	  	tree->Branch("x3_m_x4", &x3_m_x4);
	  	tree->Branch("y1_m_y2", &y1_m_y2);
	  	tree->Branch("y1_m_y3", &y1_m_y3);
	  	tree->Branch("y1_m_y4", &y1_m_y4);
	  	tree->Branch("y2_m_y3", &y2_m_y3);
	  	tree->Branch("y2_m_y4", &y2_m_y4);
	  	tree->Branch("y3_m_y4", &y3_m_y4);

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
	
	std::cout<<"Number of events for alignment with "<<nLayers<<" wire chambers: "<<acceptedCounter<<" / "<<eventCounter<<std::endl;
}

void MillepedeBinaryWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MillepedeBinaryWriter);