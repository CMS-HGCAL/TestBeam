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
	  	double time_DWC1, time_DWC2, time_DWC3, time_DWC4;
		double reco1_x, reco1_y, reco2_x, reco2_y, reco3_x, reco3_y, reco4_x, reco4_y;
		double z1, z2, z3, z4;
	  	int dwc1_z, dwc2_z, dwc3_z, dwc4_z;
	  	double x1_m_x2, x1_m_x3, x1_m_x4, x2_m_x3, x2_m_x4, x3_m_x4;
	  	double y1_m_y2, y1_m_y3, y1_m_y4, y2_m_y3, y2_m_y4, y3_m_y4;
	  	int dwc1_Ntimestamps, dwc2_Ntimestamps, dwc3_Ntimestamps, dwc4_Ntimestamps;
	  	int dwc1_goodMeasurement, dwc2_goodMeasurement, dwc3_goodMeasurement, dwc4_goodMeasurement; 
	  	int dwc1_goodMeasurementX, dwc2_goodMeasurementX, dwc3_goodMeasurementX, dwc4_goodMeasurementX; 
	  	int dwc1_goodMeasurementY, dwc2_goodMeasurementY, dwc3_goodMeasurementY, dwc4_goodMeasurementY; 
	  	int N_goodMeasurements, N_goodMeasurements_X, N_goodMeasurements_Y;
  		
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

	reco1_x = dwcs->at(0).x;
	reco1_y = dwcs->at(0).y;
	reco2_x = dwcs->at(1).x;
	reco2_y = dwcs->at(1).y;
	reco3_x = dwcs->at(2).x;
	reco3_y = dwcs->at(2).y;
	reco4_x = dwcs->at(3).x;
	reco4_y = dwcs->at(3).y;
	
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
	}

	if ((rd->event==-1)&&writeMinimal) return;
	tree->Fill();
		

}// analyze ends here


void DWC_NTupelizer::beginJob() {	
	tree = fs->make<TTree>("dwc_reco", "dwc_reco");

	tree->Branch("run", &run);
	tree->Branch("runType", &runType);
	tree->Branch("event", &n_event);
	tree->Branch("goodDWC_Measurement", &goodDWC_Measurement);
	tree->Branch("reco1_x", &reco1_x);
	tree->Branch("reco1_y", &reco1_y);
	tree->Branch("z1", &z1);
	tree->Branch("reco2_x", &reco2_x);
	tree->Branch("reco2_y", &reco2_y);
	tree->Branch("z2", &z2);
	tree->Branch("reco3_x", &reco3_x);
	tree->Branch("reco3_y", &reco3_y);
	tree->Branch("z3", &z3);
	tree->Branch("reco4_x", &reco4_x);
	tree->Branch("reco4_y", &reco4_y);
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
		tree->Branch("dwc1_goodMeasurementX", &dwc1_goodMeasurementX);
		tree->Branch("dwc2_goodMeasurementX", &dwc2_goodMeasurementX);
		tree->Branch("dwc3_goodMeasurementX", &dwc3_goodMeasurementX);
		tree->Branch("dwc4_goodMeasurementX", &dwc4_goodMeasurementX); 
		tree->Branch("dwc1_goodMeasurementY", &dwc1_goodMeasurementY);
		tree->Branch("dwc2_goodMeasurementY", &dwc2_goodMeasurementY);
		tree->Branch("dwc3_goodMeasurementY", &dwc3_goodMeasurementY);
		tree->Branch("dwc4_goodMeasurementY", &dwc4_goodMeasurementY); 
		
		tree->Branch("N_goodMeasurements", &N_goodMeasurements);
		tree->Branch("N_goodMeasurements_X", &N_goodMeasurements_X);
		tree->Branch("N_goodMeasurements_Y", &N_goodMeasurements_Y);
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