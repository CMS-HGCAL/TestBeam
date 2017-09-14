/* 
 * Determination of the position resolution of the setup.
 */

/**
	@Author: Thorben Quast <tquast>
		14 September 2017
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"  
#include <sstream>
#include <iomanip>
#include <set>
 
//#define DEBUG

class MIPFinder : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit MIPFinder(const edm::ParameterSet&);
		~MIPFinder();
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

		std::string pathToMIPWindowFile;

		std::map<int,TH2F*> m_h_rechitEnergyPerDWCE;
		std::map<int,TH1F*> m_h_rechitEnergy;
		TH2F* h_DWCE_occupancy;

		int n_bins_DWCE;
		double max_dim_x_DWCE;
		double max_dim_y_DWCE;
};

MIPFinder::MIPFinder(const edm::ParameterSet& iConfig) {	
	
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	MWCToken= consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));

	//read the configuration
	pathToMIPWindowFile = iConfig.getParameter<std::string>("pathToMIPWindowFile");
	n_bins_DWCE = iConfig.getParameter<int>("n_bins_DWCE");
	max_dim_x_DWCE = iConfig.getParameter<double>("max_dim_x_DWCE");
	max_dim_y_DWCE = iConfig.getParameter<double>("max_dim_y_DWCE");

	h_DWCE_occupancy = fs->make<TH2F>("DWC_E_Occupancy", "DWC_E_Occupancy", n_bins_DWCE, -max_dim_x_DWCE, max_dim_x_DWCE, n_bins_DWCE, -max_dim_y_DWCE, max_dim_y_DWCE);		//DWC dimension and binning to be configured

	TH2F* htmp2;
	TH1F* htmp1;
	std::ostringstream os( std::ostringstream::ate );
	for(size_t ib = 1; ib<=HGCAL_TB_GEOMETRY::NUMBER_OF_HEXABOARD; ib++) {
		for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
			os.str("");
			os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
			TFileDirectory chip_dir = fs->mkdir( os.str().c_str() );
			for( size_t ichan=0; ichan<=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
				if (ichan%2==1) continue;
				os.str("");
				os << "Channel" << ichan;
				TFileDirectory channel_dir = chip_dir.mkdir( os.str().c_str() );
				
				os.str("");
				os << "EnergyVsDWCE_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp2=channel_dir.make<TH2F>("RechitEnergyVsDWCE", os.str().c_str(), n_bins_DWCE, -max_dim_x_DWCE, max_dim_x_DWCE, n_bins_DWCE, -max_dim_y_DWCE, max_dim_y_DWCE);
				m_h_rechitEnergyPerDWCE.insert( std::pair<int,TH2F*>(ib*1000+iski*100+ichan, htmp2) );
				
				os.str("");
				os << "Energy_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp1=channel_dir.make<TH1F>("RechitEnergy",os.str().c_str(), 149, 1., 150.);
				m_h_rechitEnergy.insert( std::pair<int,TH1F*>(ib*1000+iski*100+ichan, htmp1) );
				
			}
		}
	}


}

MIPFinder::~MIPFinder() {
	return;
}

// ------------ method called for each event  ------------
void MIPFinder::analyze(const edm::Event& event, const edm::EventSetup& setup) {

	edm::Handle<RunData> rd;
 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	
	int evId = event.id().event();
	int run = rd->run;
	int pdgID = std::atoi( (rd->runType).c_str() );
	double energy = rd->energy;
	
	if (rd->hasDanger) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because somthing went wrong"<<std::endl;
		return;
	}

	if (pdgID != 13) {
		std::cout<<"Run is not a dedicated muon run."<<std::endl;
		return;
	}
	#ifdef DEBUG
		int eventCounter = rd->event;
		std::cout<<"run: "<<run<<"  energy: "<<energy<<"  pdgID:" << pdgID<<"   eventCounter: "<<eventCounter<<std::endl;
	#endif


	//Obtain the wire chamber information
	edm::Handle<WireChambers> dwcs;
	event.getByToken(MWCToken, dwcs);
		
	if (!dwcs->at(0).goodMeasurement) {
		return;
	}
	double DWCE_x = dwcs->at(0).x;
	double DWCE_y = dwcs->at(0).y;
	#ifdef DEBUG
		std::cout<<dwcs->at(0).x<<"  "<<dwcs->at(0).y<<"   "<<dwcs->at(0).z<<"  "<<dwcs->at(0).goodMeasurement<<std::endl;
	#endif
	h_DWCE_occupancy->Fill(DWCE_x, DWCE_y);

	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);

	
	//fill the rechits:
	for(auto Rechit : *Rechits) {	
		int layer = (Rechit.id()).layer();
		int skiroc = Rechit.skiroc() % 4;
		int channel = Rechit.channel();
  		double energyHG = Rechit.energyHigh();
  		int key = layer*1000+skiroc*100+channel;
  		
  		m_h_rechitEnergyPerDWCE[key]->Fill(DWCE_x, DWCE_y, energyHG);
  		m_h_rechitEnergy[key]->Fill(energyHG);		
	}
	
}// analyze ends here

void MIPFinder::beginJob() {	
}

void MIPFinder::endJob() {
	
}

void MIPFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MIPFinder);