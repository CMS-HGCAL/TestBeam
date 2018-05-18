/* 
 * Computes reconstructed energy spectra as a function of the impact position in the coordinate frame of DATURA telescope tracks
 */

/**
	@Author: Thorben Quast <tquast>
		18th May 2018
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
#include "HGCal/DataFormats/interface/HGCalTBDATURATelescopeData.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
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

class DATURATelescopeCorrelator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit DATURATelescopeCorrelator(const edm::ParameterSet&);
		~DATURATelescopeCorrelator();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;
		

		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;	 	
		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<std::vector<HGCalTBDATURATelescopeData> > DATURATrackToken;		

		std::string m_electronicMap;
		std::string m_detectorLayoutFile;
		struct {
			HGCalElectronicsMap emap_;
			HGCalTBDetectorLayout layout_;
		} essource_;

		int m_NHexaBoards;

		int n_bins;
		double max_dim_x_DUT;
		double max_dim_y_DUT;		
	  	
		//output histograms
		std::map<int,TH2F*> m_h_board_occupancy;		
		std::map<int,TH2F*> m_h_rechitEnergyPerDUTAveraged;
		std::map<int,TH2F*> m_h_rechitEfficiencyPerDUT;
		
};

DATURATelescopeCorrelator::DATURATelescopeCorrelator(const edm::ParameterSet& iConfig) {	
	
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	DATURATrackToken= consumes<std::vector<HGCalTBDATURATelescopeData> >(iConfig.getParameter<edm::InputTag>("DATURATelescopeData"));
	
	m_electronicMap = iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt");
	m_detectorLayoutFile = iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt");
	m_NHexaBoards= iConfig.getUntrackedParameter<int>("NHexaBoards", 10);


	n_bins = iConfig.getParameter<int>("n_bins");
	max_dim_x_DUT = iConfig.getParameter<double>("max_dim_x_DUT");
	max_dim_y_DUT = iConfig.getParameter<double>("max_dim_y_DUT");




	

	TH2F* htmp2;
	
	std::ostringstream os( std::ostringstream::ate );
	for(int ib = 0; ib<m_NHexaBoards; ib++) {
		os.str("");
		os << "Occupancy_Board"<<ib;
		htmp2=fs->make<TH2F>(os.str().c_str(), os.str().c_str(), n_bins, -max_dim_x_DUT, max_dim_x_DUT, n_bins, -max_dim_y_DUT, max_dim_y_DUT);
		m_h_board_occupancy.insert( std::pair<int,TH2F*>(ib, htmp2) );
	
		for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
			os.str("");
			os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
			TFileDirectory chip_dir = fs->mkdir( os.str().c_str() );
			for( size_t ichan=0; ichan<=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
				if (ichan%2==1) continue;
				int key = ib*1000+iski*100+ichan;

				os.str("");
				os << "Channel" << ichan;
				TFileDirectory channel_dir = chip_dir.mkdir( os.str().c_str() );
				
		
				os.str("");
				os << "AverageRechitEnergyPerDUT"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp2 = channel_dir.make<TH2F>("AverageRechitEnergyPerDUT", os.str().c_str(), n_bins, -max_dim_x_DUT, max_dim_x_DUT, n_bins, -max_dim_y_DUT, max_dim_y_DUT);
				m_h_rechitEnergyPerDUTAveraged.insert( std::pair<int,TH2F*>(key, htmp2));

	
				os.str("");
				os << "HGRechitEfficiencyPerDUT"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp2=channel_dir.make<TH2F>("HGRechitEfficiency", os.str().c_str(), n_bins, -max_dim_x_DUT, max_dim_x_DUT, n_bins, -max_dim_y_DUT, max_dim_y_DUT);
				m_h_rechitEfficiencyPerDUT.insert( std::pair<int,TH2F*>(key, htmp2) );	

			}
		}
	}

	


}

DATURATelescopeCorrelator::~DATURATelescopeCorrelator() {
	return;
}

// ------------ method called for each event  ------------
void DATURATelescopeCorrelator::analyze(const edm::Event& event, const edm::EventSetup& setup) {

	edm::Handle<RunData> rd;
 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	
	int evId = event.id().event();
	int run = rd->run;
	double energy = rd->energy;
	
	if (rd->booleanUserRecords.has("hasDanger")&&rd->booleanUserRecords.get("hasDanger")) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because somthing went wrong"<<std::endl;
		return;
	}

	if (!rd->booleanUserRecords.has("hasValidDATURAMeasurement")) {
		std::cout<<"No DATURA data tagged to the run data "<<std::endl;
		return;
	}
	if (!rd->booleanUserRecords.get("hasValidDATURAMeasurement")) return;

	//load the DATURA tracks
	edm::Handle<std::vector<HGCalTBDATURATelescopeData> > daturatracks;	
	event.getByToken(DATURATrackToken, daturatracks);

	//load the rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);	

	int nTrackCounter=0;
	for(auto daturatrack : *daturatracks) {
		nTrackCounter++;
		//maximum two tracks per event
		//chi2 based selection

		for(int ib = 0; ib<m_NHexaBoards; ib++) {
			int layer=essource_.layout_.getLayerWithModuleIndex(ib).layerID()+1;
			if (daturatrack.Extrapolation_XY_Chi2(layer).first > 100.) continue;
			if (daturatrack.Extrapolation_XY_Chi2(layer).second > 100.) continue;
			m_h_board_occupancy[ib]->Fill(daturatrack.Extrapolation_XY(layer).first, daturatrack.Extrapolation_XY(layer).second);

			for(auto Rechit : *Rechits) {	
				int rechit_layer = (Rechit.id()).layer();
				if (rechit_layer!=layer) continue;

		  		float DUT_x = daturatrack.Extrapolation_XY(layer).first;
		  		float DUT_y = daturatrack.Extrapolation_XY(layer).second;
				
				HGCalTBElectronicsId eid( essource_.emap_.detId2eid( Rechit.id().rawId() ) );
				int skiroc = eid.iskiroc_rawhit();
				int board = skiroc / 4;
				
				skiroc = eid.iskiroc_rawhit() % 4;
				int channel = eid.ichan();
		  		int key = board*1000+skiroc*100+channel;
		  		double energy = (Rechit.checkFlag(HGCalTBRecHit::kGood)&&(!Rechit.checkFlag(HGCalTBRecHit::kHighGainSaturated))) ? Rechit.energyHigh() : 0;

		  		m_h_rechitEnergyPerDUTAveraged[key]->Fill(DUT_x, DUT_y, energy);
				int hit = (energy > 0) ? 1: 0;
				m_h_rechitEfficiencyPerDUT[key]->Fill(DUT_x, DUT_y, hit);
			}
		}
	}
}

void DATURATelescopeCorrelator::beginJob() {	
	
	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(m_electronicMap);
	if (!io.load(fip.fullPath(), essource_.emap_)) {
		throw cms::Exception("Unable to load electronics map");
	};

  	fip=edm::FileInPath(m_detectorLayoutFile);
	if (!io.load(fip.fullPath(), essource_.layout_)) {
		throw cms::Exception("Unable to load detector layout file");
	};
}

void DATURATelescopeCorrelator::endJob() {
	for(int ib = 0; ib<m_NHexaBoards; ib++) {
		for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
			for( size_t ichan=0; ichan<=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
				if (ichan%2==1) continue;
				int key = ib*1000+iski*100+ichan;
				m_h_rechitEnergyPerDUTAveraged[key]->Divide(m_h_board_occupancy[ib]);
				m_h_rechitEnergyPerDUTAveraged[key]->SetStats(false);
				m_h_rechitEnergyPerDUTAveraged[key]->GetZaxis()->SetRangeUser(0., 120.);		//skiroc2-cms: usual scale for MIPs are ~50 HG ADC  	

				m_h_rechitEfficiencyPerDUT[key]->Divide(m_h_board_occupancy[ib]);
				m_h_rechitEfficiencyPerDUT[key]->SetStats(false);
				m_h_rechitEfficiencyPerDUT[key]->GetZaxis()->SetRangeUser(0., 1.);		//skiroc2-cms: usual scale for MIPs are ~50 HG ADC 
			}
		}
	}	
}

void DATURATelescopeCorrelator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(DATURATelescopeCorrelator);