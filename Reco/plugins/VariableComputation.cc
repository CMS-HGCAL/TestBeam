/* 
 * Computation of variables.
 */

/**
	@Author: Thorben Quast <tquast>
		22 November 2017
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
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBCommonModeNoise.h"
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"  
#include <sstream>
#include <fstream>
#include <iomanip>
#include <set>


//#define DEBUG

class VariableComputation : public edm::EDProducer {
	public:
		explicit VariableComputation(const edm::ParameterSet&);
		~VariableComputation();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob();
		virtual void produce(edm::Event& , const edm::EventSetup&);
		virtual void endJob();


		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;	 	
		edm::EDGetTokenT<std::map<int, commonModeNoise> > CommonModeNoiseMap_Token;	 	
		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<WireChambers> DWCToken;		
		edm::EDGetTokenT<HGCalTBDWCTrack> DWCTrackToken;		

		std::string m_UserRecordCollectionName;

		std::string m_electronicMap;
		std::string m_detectorLayoutFile;
		struct {
			HGCalElectronicsMap emap_;
			HGCalTBDetectorLayout layout_;
		} essource_;

		int m_NHexaBoards;
};

VariableComputation::VariableComputation(const edm::ParameterSet& iConfig) {	
	
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	CommonModeNoiseMap_Token = consumes<std::map<int, commonModeNoise>>(iConfig.getParameter<edm::InputTag>("HGCALTBCOMMONMODENOISE"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	DWCToken= consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
	DWCTrackToken= consumes<HGCalTBDWCTrack>(iConfig.getParameter<edm::InputTag>("DWCTRACKS"));
	
	m_UserRecordCollectionName = iConfig.getUntrackedParameter<std::string>("UserRecordCollectionName","DoubleUserRecords");
	
	m_electronicMap = iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt");
	m_detectorLayoutFile = iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt");
	m_NHexaBoards= iConfig.getUntrackedParameter<int>("NHexaBoards", 10);


	produces <UserRecords<double> >(m_UserRecordCollectionName);
}

VariableComputation::~VariableComputation() {
	return;
}

// ------------ method called for each event  ------------
void VariableComputation::produce(edm::Event& event, const edm::EventSetup& setup) {
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);
	
	edm::Handle<std::map<int, commonModeNoise>> cmMap;
	//event.getByToken(CommonModeNoiseMap_Token, cmMap);
	
	edm::Handle<HGCalTBDWCTrack> dwctrack;
	event.getByToken(DWCTrackToken, dwctrack);

	edm::Handle<WireChambers> dwcs;
	event.getByToken(DWCToken, dwcs);

	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);

	std::auto_ptr<UserRecords<double> > UR(new UserRecords<double>);
	
	UR->add("eventID", rd->event);
	UR->add("run", rd->run);
	UR->add("pdgID", rd->pdgID);
	UR->add("beamEnergy", rd->energy);
	UR->add("configuration", rd->configuration);
	UR->add("runType", rd->runType);

	event.put(UR, m_UserRecordCollectionName);
	
}// analyze ends here

void VariableComputation::beginJob() {	
}

void VariableComputation::endJob() {
}

void VariableComputation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VariableComputation);