/* 
 * Write out residuals and derivatives to pass forward to millepede. Relies on the fact that all layers have proper hits. 
 */

/**
	@Author: Thorben Quast <tquast>
		09 Nov 2016
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
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/Reco/interface/cnpy.h"

class NumpyConverter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit NumpyConverter(const edm::ParameterSet&);
		~NumpyConverter();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------

		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;	 	
		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<std::map<int, WireChamberData> > DWCToken;		
		edm::EDGetTokenT<HGCalTBDWCTrack> DWCTrackToken;		

		std::string m_electronicMap;
		struct {
			HGCalElectronicsMap emap_;
		} essource_;

		//relevant for rechits storage
		uint m_NHexaBoards;
    	uint m_Sensorsize;
    	
    	//relevant for DWC storage
    	uint m_NLayers;
    	uint m_NDWCs;

    	std::string m_outputFilePath;
    	std::string m_eventDataIdentifier;
    	std::string m_rechitIdentifier;
    	std::string m_dwcIdentifier;
    	std::string m_dwcReferenceIdentifier;

    	//delimiters of the hexagonal coordinate system
    	int u_max;
    	int u_min;
    	uint range_u;
    	int v_max;
    	int v_min;
    	uint range_v;


    	uint Nevents;
    	std::vector<float> event_data;
    	std::vector<float> rechit_data;
    	std::vector<float> dwc_data;
    	std::vector<float> dwc_track_data;
};

NumpyConverter::NumpyConverter(const edm::ParameterSet& iConfig) {
	usesResource("TFileService");

	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	DWCToken= consumes<std::map<int, WireChamberData> >(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
	DWCTrackToken= consumes<HGCalTBDWCTrack>(iConfig.getParameter<edm::InputTag>("DWCTRACKS"));
	
	m_electronicMap = iConfig.getUntrackedParameter<std::string>("electronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt");
	m_NHexaBoards= iConfig.getUntrackedParameter<uint>("NHexaBoards", 10);
    m_Sensorsize=iConfig.getUntrackedParameter<uint>("Sensorsize", 128);
    
    m_NLayers=iConfig.getUntrackedParameter<uint>("NLayers", 10);
    m_NDWCs=iConfig.getUntrackedParameter<uint>("NDWCs", 4);

    m_outputFilePath=iConfig.getParameter<std::string>("outputFilePath");
	m_eventDataIdentifier=iConfig.getUntrackedParameter<std::string>("eventDataIdentifier", "event");
	m_rechitIdentifier=iConfig.getUntrackedParameter<std::string>("rechitIdentifier", "rechits");
	m_dwcIdentifier=iConfig.getUntrackedParameter<std::string>("dwcIdentifier", "dwcs");
	m_dwcReferenceIdentifier=iConfig.getUntrackedParameter<std::string>("dwcReferenceIdentifier", "dwcReference");

    Nevents = 0;

    if (m_Sensorsize==128) {
    	u_max = 7;
    	u_min = -7;
    	v_max = 7;
    	v_min = -7;
    } else {		//other geometries to be implemented
      	u_max = 7;
    	u_min = -7;
    	v_max = 7;
    	v_min = -7;  	
    }
    range_u = u_max-u_min + 1;
    range_v = v_max-v_min + 1;
}//constructor ends here

NumpyConverter::~NumpyConverter() {
	return;
}

// ------------ method called for each event  ------------
void NumpyConverter::analyze(const edm::Event& event, const edm::EventSetup& setup) {
	Nevents++;
	


	//fill the event info
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);
	event_data.push_back((float)event.id().event());
	event_data.push_back((float)rd->run);
	event_data.push_back((float)rd->pdgID);
	event_data.push_back((float)rd->energy);
	
	
 	//fill the rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);
	//initialize data with zeros
	float**** data = new float***[m_NHexaBoards];
	for (uint b=0; b<m_NHexaBoards; b++) {
		data[b] = new float**[range_u];
		for (uint u=0; u<(range_u); u++) {
			data[b][u] = new float*[range_v];
			for (uint v=0; v<(range_v); v++) {
				data[b][u][v] = new float[2];
				data[b][u][v][0] = data[b][u][v][1] = 0.;
			}
		}
	}
	for(auto Rechit : *Rechits) {	
		HGCalTBElectronicsId eid( essource_.emap_.detId2eid( Rechit.id().rawId() ) );
		HGCalTBDetId detId = HGCalTBDetId(Rechit.id().rawId());
		int b = eid.iskiroc_rawhit() / 4;
		int u = detId.iu()-u_min;
		int v = detId.iv()-v_min;
		data[b][u][v][0] = Rechit.energy();
		data[b][u][v][1] = Rechit.time();
	}
	for (uint b=0; b<m_NHexaBoards; b++) for (uint u=0; u<(range_u); u++) for (uint v=0; v<(range_v); v++) {		
		rechit_data.push_back(data[b][u][v][0]);
		rechit_data.push_back(data[b][u][v][1]);
	}
	for (uint b=0; b<m_NHexaBoards; b++) {
		for (uint u=0; u<(range_u); u++) {
			for (uint v=0; v<(range_v); v++) {
				delete[] data[b][u][v];
			}
			delete[] data[b][u];
		}
		delete[] data[b];
	}		
	delete[] data;	
	

	//fill the DWC information
	edm::Handle<std::map<int, WireChamberData> > dwcs;
	event.getByToken(DWCToken, dwcs);
	for (int d=0; d<(int)m_NDWCs; d++) {
		dwc_data.push_back((float)dwcs->at(d).x);
		dwc_data.push_back((float)dwcs->at(d).y);
		dwc_data.push_back((float)dwcs->at(d).z);
	}

	//fill the DWC track information
	edm::Handle<HGCalTBDWCTrack> dwctrack;
	event.getByToken(DWCTrackToken, dwctrack);
	if (dwctrack->valid){
		dwc_track_data.push_back((float)dwctrack->referenceType);
		dwc_track_data.push_back((float)dwctrack->referenceType);
		for (int l=1; l<=(int)m_NLayers; l++) {
			dwc_track_data.push_back((float)dwctrack->DWCExtrapolation_XY(l).first);
			dwc_track_data.push_back((float)dwctrack->DWCExtrapolation_XY(l).second);
		}
		dwc_track_data.push_back((float)dwctrack->chi2_x);
		dwc_track_data.push_back((float)dwctrack->chi2_y);
	} else {
		for (uint i=0; i<2*(2+m_NLayers); i++) dwc_track_data.push_back(-1.);
	}

}// analyze ends here


void NumpyConverter::beginJob() {	
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(m_electronicMap);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	};
}

void NumpyConverter::endJob() {
	cnpy::npz_save(m_outputFilePath, m_eventDataIdentifier.c_str(), &event_data[0],{Nevents, 4}, "w");
	cnpy::npz_save(m_outputFilePath, m_rechitIdentifier.c_str(), &rechit_data[0],{Nevents, m_NHexaBoards, range_u, range_v, 2}, "a");
	cnpy::npz_save(m_outputFilePath, m_dwcIdentifier.c_str(), &dwc_data[0],{Nevents, m_NDWCs, 3}, "a");
	cnpy::npz_save(m_outputFilePath, m_dwcReferenceIdentifier.c_str(), &dwc_track_data[0],{Nevents, 1+m_NLayers+1, 2}, "a");
}

void NumpyConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NumpyConverter);