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
#include "HGCal/Reco/interface/Tracks.h"
#include "HGCal/Reco/interface/Sensors.h"

#include "TTree.h"
#include "TFile.h"


class DWC_NTupelizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit DWC_NTupelizer(const edm::ParameterSet&);
		~DWC_NTupelizer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------

		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<WireChambers> MWCToken;
		edm::Service<TFileService> fs;
		bool writeMinimal;

  	int run, runType, n_event, goodDWC_Measurement;
  	double triggerTimeDiff;
  	double time_DWC1, time_DWC2, time_DWC3, time_DWC4;
		double reco1_x, reco1_y, reco2_x, reco2_y, reco3_x, reco3_y, reco4_x, reco4_y;
		double res1_x, res1_y, res2_x, res2_y, res3_x, res3_y, res4_x, res4_y;
		double z1, z2, z3, z4;
  	int dwc1_z, dwc2_z, dwc3_z, dwc4_z;
  	double x1_m_x2, x1_m_x3, x1_m_x4, x2_m_x3, x2_m_x4, x3_m_x4;
  	double y1_m_y2, y1_m_y3, y1_m_y4, y2_m_y3, y2_m_y4, y3_m_y4;
  	int dwc1_Ntimestamps, dwc2_Ntimestamps, dwc3_Ntimestamps, dwc4_Ntimestamps;
  	int dwc1_goodMeasurement, dwc2_goodMeasurement, dwc3_goodMeasurement, dwc4_goodMeasurement; 
  	int dwc1_goodMeasurementX, dwc2_goodMeasurementX, dwc3_goodMeasurementX, dwc4_goodMeasurementX; 
  	int dwc1_goodMeasurementY, dwc2_goodMeasurementY, dwc3_goodMeasurementY, dwc4_goodMeasurementY; 
  	int N_goodMeasurements, N_goodMeasurements_X, N_goodMeasurements_Y;
  	double residual1_x, residual1_y, residual2_x, residual2_y, residual3_x, residual3_y, residual4_x, residual4_y; 


  	int NDF_x; double chi2_x, referenceError_x;
  	int NDF_y; double chi2_y, referenceError_y;
		std::map<int, SensorHitMap*> Sensors;
		ParticleTrack* Track;  		
		TTree* tree;

};

DWC_NTupelizer::DWC_NTupelizer(const edm::ParameterSet& iConfig) {
	usesResource("TFileService");

	// initialization	
	MWCToken = consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
	RunDataToken = consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));

	writeMinimal = iConfig.getParameter<bool>("writeMinimal");

	tree = NULL;


}//constructor ends here

DWC_NTupelizer::~DWC_NTupelizer() {
	return;
}

// ------------ method called for each event  ------------
void DWC_NTupelizer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
 	//get the relevant event information
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);

	//get the multi wire chambers
	edm::Handle<WireChambers> dwcs;
	event.getByToken(MWCToken, dwcs);
	
	run = rd->run;
	runType = std::atoi((rd->runType).c_str());
	n_event = rd->event;
	goodDWC_Measurement = rd->hasValidMWCMeasurement ? 1 : 0;
	triggerTimeDiff = rd->triggerDeltaT_to_TDC;

	reco1_x = dwcs->at(0).x;
	reco1_y = dwcs->at(0).y;
	reco2_x = dwcs->at(1).x;
	reco2_y = dwcs->at(1).y;
	reco3_x = dwcs->at(2).x;
	reco3_y = dwcs->at(2).y;
	reco4_x = dwcs->at(3).x;
	reco4_y = dwcs->at(3).y;

	res1_x = dwcs->at(0).res_x;
	res1_y = dwcs->at(0).res_y;
	res2_x = dwcs->at(1).res_x;
	res2_y = dwcs->at(1).res_y;
	res3_x = dwcs->at(2).res_x;
	res3_y = dwcs->at(2).res_y;
	res4_x = dwcs->at(3).res_x;
	res4_y = dwcs->at(3).res_y;
	
	z1 = dwcs->at(0).z;
	z2 = dwcs->at(1).z;
	z3 = dwcs->at(2).z;
	z4 = dwcs->at(3).z;

	dwc1_goodMeasurement = dwcs->at(0).goodMeasurement;
	dwc2_goodMeasurement = dwcs->at(1).goodMeasurement;
	dwc3_goodMeasurement = dwcs->at(2).goodMeasurement;
	dwc4_goodMeasurement = dwcs->at(3).goodMeasurement;

	if (!writeMinimal) {	
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


		dwc1_Ntimestamps = dwcs->at(0).recordedTimeStamps;
		dwc2_Ntimestamps = dwcs->at(1).recordedTimeStamps;
		dwc3_Ntimestamps = dwcs->at(2).recordedTimeStamps;
		dwc4_Ntimestamps = dwcs->at(3).recordedTimeStamps;
		
		dwc1_goodMeasurementX = dwcs->at(0).goodMeasurement_X;
		dwc2_goodMeasurementX = dwcs->at(1).goodMeasurement_X;
		dwc3_goodMeasurementX = dwcs->at(2).goodMeasurement_X;
		dwc4_goodMeasurementX = dwcs->at(3).goodMeasurement_X;
		
		dwc1_goodMeasurementY = dwcs->at(0).goodMeasurement_Y;
		dwc2_goodMeasurementY = dwcs->at(1).goodMeasurement_Y;
		dwc3_goodMeasurementY = dwcs->at(2).goodMeasurement_Y;
		dwc4_goodMeasurementY = dwcs->at(3).goodMeasurement_Y;

		N_goodMeasurements = N_goodMeasurements_X = N_goodMeasurements_Y = 0;

		if(dwc1_goodMeasurement) N_goodMeasurements++;
		if(dwc2_goodMeasurement) N_goodMeasurements++;
		if(dwc3_goodMeasurement) N_goodMeasurements++;
		if(dwc4_goodMeasurement) N_goodMeasurements++;

		if(dwc1_goodMeasurementX) N_goodMeasurements_X++;
		if(dwc2_goodMeasurementX) N_goodMeasurements_X++;
		if(dwc3_goodMeasurementX) N_goodMeasurements_X++;
		if(dwc4_goodMeasurementX) N_goodMeasurements_X++;
		
		if(dwc1_goodMeasurementY) N_goodMeasurements_Y++;
		if(dwc2_goodMeasurementY) N_goodMeasurements_Y++;
		if(dwc3_goodMeasurementY) N_goodMeasurements_Y++;
		if(dwc4_goodMeasurementY) N_goodMeasurements_Y++;



		//fit for x
		if (N_goodMeasurements_X>=2) {
			//first layer of HGCal as reference
			Track = new ParticleTrack();
			Sensors[5] = new SensorHitMap(5);
			Sensors[5]->setLabZ(0., 0);
			Sensors[5]->setParticleEnergy(rd->energy);
			Track->addReferenceSensor(Sensors[5]);
			
			for (size_t n_layer=0; n_layer<4; n_layer++) {	
				Sensors[n_layer] = new SensorHitMap(n_layer);				
				Sensors[n_layer]->setLabZ(dwcs->at(n_layer).z, 0.);
				Sensors[n_layer]->setCenterHitPosition(dwcs->at(n_layer).x, dwcs->at(n_layer).y , dwcs->at(n_layer).res_x, dwcs->at(n_layer).res_y);
				//std::cout<<n_layer<<": "<<dwcs->at(n_layer).y<<" good measurement ? "<<dwcs->at(n_layer).goodMeasurement_Y<<std::endl;

				Sensors[n_layer]->setParticleEnergy(rd->energy);
				//Sensors[n_layer]->setResidualResolution(dwcs->at(n_layer).res_x);	
				Track->addFitPoint(Sensors[n_layer]);
			}
			
			Track->fitTrack(LINEFITANALYTICAL);
			NDF_x = Track->getNDF(1);
			chi2_x = Track->getChi2(1);
			referenceError_x = Track->calculateReferenceErrorXY().first;
			residual1_x = Sensors[0]->getHitPosition().first-Track->calculatePositionXY(Sensors[0]->getLabZ(), 0).first;
			residual2_x = Sensors[1]->getHitPosition().first-Track->calculatePositionXY(Sensors[1]->getLabZ(), 0).first;
			residual3_x = Sensors[2]->getHitPosition().first-Track->calculatePositionXY(Sensors[2]->getLabZ(), 0).first;
			residual4_x = Sensors[3]->getHitPosition().first-Track->calculatePositionXY(Sensors[3]->getLabZ(), 0).first;

			for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
				delete (*it).second;
			};	Sensors.clear();
			delete Track;
		} else {
			NDF_x = 0;
			chi2_x = -999;
			referenceError_x = -999;
			residual1_x = residual2_x = residual3_x = residual4_x = -999;
		}

		//fit for y
		if (N_goodMeasurements_Y>=2) {
			//first layer of HGCal as reference
			Track = new ParticleTrack();
			Sensors[5] = new SensorHitMap(5);
			Sensors[5]->setLabZ(0., 0);
			Sensors[5]->setParticleEnergy(rd->energy);
			Track->addReferenceSensor(Sensors[5]);
			
			for (size_t n_layer=0; n_layer<4; n_layer++) {	
				Sensors[n_layer] = new SensorHitMap(n_layer);				
				Sensors[n_layer]->setLabZ(dwcs->at(n_layer).z, 0.);
				Sensors[n_layer]->setCenterHitPosition(dwcs->at(n_layer).x, dwcs->at(n_layer).y , dwcs->at(n_layer).res_x, dwcs->at(n_layer).res_y);
				//std::cout<<n_layer<<": "<<dwcs->at(n_layer).y<<" good measurement ? "<<dwcs->at(n_layer).goodMeasurement_Y<<std::endl;

				Sensors[n_layer]->setParticleEnergy(rd->energy);
				//Sensors[n_layer]->setResidualResolution(dwcs->at(n_layer).res_x);	
				Track->addFitPoint(Sensors[n_layer]);
			}
			
			Track->fitTrack(LINEFITANALYTICAL);
			NDF_y = Track->getNDF(2);
			chi2_y = Track->getChi2(2);
			referenceError_y = Track->calculateReferenceErrorXY().second;
			residual1_y = Sensors[0]->getHitPosition().second-Track->calculatePositionXY(Sensors[0]->getLabZ(), 0).second;
			residual2_y = Sensors[1]->getHitPosition().second-Track->calculatePositionXY(Sensors[1]->getLabZ(), 0).second;
			residual3_y = Sensors[2]->getHitPosition().second-Track->calculatePositionXY(Sensors[2]->getLabZ(), 0).second;
			residual4_y = Sensors[3]->getHitPosition().second-Track->calculatePositionXY(Sensors[3]->getLabZ(), 0).second;


			for (std::map<int, SensorHitMap*>::iterator it=Sensors.begin(); it!=Sensors.end(); it++) {
				delete (*it).second;
			};	Sensors.clear();
			delete Track;

		} else {
			NDF_y = 0;
			chi2_y = -999;
			referenceError_y = -999;
			residual1_y = residual2_y = residual3_y = residual4_y = -999;
		}


	}

	if ((rd->event==-1)&&writeMinimal) return;
	tree->Fill();
		

}// analyze ends here


void DWC_NTupelizer::beginJob() {	
	tree = fs->make<TTree>("dwc_reco", "dwc_reco");

	tree->Branch("run", &run);
	tree->Branch("runType", &runType);
	tree->Branch("event", &n_event);
	tree->Branch("triggerTimeDifference", &triggerTimeDiff);
	tree->Branch("goodDWC_Measurement", &goodDWC_Measurement);
	tree->Branch("reco1_x", &reco1_x);
	tree->Branch("reco1_y", &reco1_y);
	tree->Branch("res1_x", &res1_x);
	tree->Branch("res1_y", &res1_y);
	tree->Branch("z1", &z1);
	tree->Branch("reco2_x", &reco2_x);
	tree->Branch("reco2_y", &reco2_y);
	tree->Branch("res2_x", &res2_x);
	tree->Branch("res2_y", &res2_y);
	tree->Branch("z2", &z2);
	tree->Branch("reco3_x", &reco3_x);
	tree->Branch("reco3_y", &reco3_y);
	tree->Branch("res3_x", &res3_x);
	tree->Branch("res3_y", &res3_y);
	tree->Branch("z3", &z3);
	tree->Branch("reco4_x", &reco4_x);
	tree->Branch("reco4_y", &reco4_y);
	tree->Branch("res4_x", &res4_x);
	tree->Branch("res4_y", &res4_y);
	tree->Branch("z4", &z4);
	tree->Branch("dwc1_goodMeasurement", &dwc1_goodMeasurement);
	tree->Branch("dwc2_goodMeasurement", &dwc2_goodMeasurement);
	tree->Branch("dwc3_goodMeasurement", &dwc3_goodMeasurement);
	tree->Branch("dwc4_goodMeasurement", &dwc4_goodMeasurement); 


	if (!writeMinimal) {
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

		tree->Branch("dwc1_Ntimestamps", &dwc1_Ntimestamps);
		tree->Branch("dwc2_Ntimestamps", &dwc2_Ntimestamps);
		tree->Branch("dwc3_Ntimestamps", &dwc3_Ntimestamps);
		tree->Branch("dwc4_Ntimestamps", &dwc4_Ntimestamps);
		tree->Branch("dwc1_goodMeasurement_x", &dwc1_goodMeasurementX);
		tree->Branch("dwc2_goodMeasurement_x", &dwc2_goodMeasurementX);
		tree->Branch("dwc3_goodMeasurement_x", &dwc3_goodMeasurementX);
		tree->Branch("dwc4_goodMeasurement_x", &dwc4_goodMeasurementX); 
		tree->Branch("dwc1_goodMeasurement_y", &dwc1_goodMeasurementY);
		tree->Branch("dwc2_goodMeasurement_y", &dwc2_goodMeasurementY);
		tree->Branch("dwc3_goodMeasurement_y", &dwc3_goodMeasurementY);
		tree->Branch("dwc4_goodMeasurement_y", &dwc4_goodMeasurementY); 
		
		tree->Branch("N_goodMeasurements", &N_goodMeasurements);
		tree->Branch("N_goodMeasurements_x", &N_goodMeasurements_X);
		tree->Branch("N_goodMeasurements_y", &N_goodMeasurements_Y);
	
		tree->Branch("NDF_x", &NDF_x);
		tree->Branch("chi2_x", &chi2_x);
		tree->Branch("referenceError_x", &referenceError_x);

		tree->Branch("NDF_y", &NDF_y);
		tree->Branch("chi2_y", &chi2_y);
		tree->Branch("referenceError_y", &referenceError_y);

		tree->Branch("residual1_x", &residual1_x);
		tree->Branch("residual1_y", &residual1_y);
		tree->Branch("residual2_x", &residual2_x);
		tree->Branch("residual2_y", &residual2_y);
		tree->Branch("residual3_x", &residual3_x);
		tree->Branch("residual3_y", &residual3_y);
		tree->Branch("residual4_x", &residual4_x);
		tree->Branch("residual4_y", &residual4_y);
	}	
}

void DWC_NTupelizer::endJob() {

}

void DWC_NTupelizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DWC_NTupelizer);