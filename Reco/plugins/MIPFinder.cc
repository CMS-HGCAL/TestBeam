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
#include "HGCal/DataFormats/interface/HGCalTBDWCTrack.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"  
#include <sstream>
#include <fstream>
#include <iomanip>
#include <set>
 

typedef std::map<int, std::vector<double> > WindowMap;

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
		void ReadDWCWindows(int);
		void ReadCurrentDWCWindows(int);


		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;	 	
		edm::EDGetTokenT<RunData> RunDataToken;	
		edm::EDGetTokenT<WireChambers> DWCToken;		
		edm::EDGetTokenT<HGCalTBDWCTrack> DWCTrackToken;		

		int m_NHexaBoards;

		std::vector<std::string> pathsToMIPWindowFiles;
	  	std::map<std::pair<int, int> ,WindowMap  >loadedDWCWindows;
		WindowMap currentDWCWindows;


		std::map<int,TH2F*> h_DUT_occupancy;
		
		std::map<int,TH2F*> m_h_rechitEnergyPerDUT;
		std::map<int,TH1F*> m_h_rechitEnergy;
		
		std::map<int,TH2F*> m_h_rechitEnergyPerDUT_selected;
		std::map<int,TH1F*> m_h_rechitEnergy_selected;

		std::map<int,TH2F*> m_h_DUT_occupancy_forBoard;		
		std::map<int,TH2F*> m_h_rechitEfficiencyPerDUT;		



		int n_bins_DWCE;
		double max_dim_x_DUT;
		double max_dim_y_DUT;
};

MIPFinder::MIPFinder(const edm::ParameterSet& iConfig) {	
	
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	DWCToken= consumes<WireChambers>(iConfig.getParameter<edm::InputTag>("MWCHAMBERS"));
	DWCTrackToken= consumes<HGCalTBDWCTrack>(iConfig.getParameter<edm::InputTag>("DWCTRACKS"));
	m_NHexaBoards= iConfig.getUntrackedParameter<int>("NHexaBoards", 10);

	//read the configuration
	pathsToMIPWindowFiles = iConfig.getParameter<std::vector<std::string> >("pathsToMIPWindowFiles");
	n_bins_DWCE = iConfig.getParameter<int>("n_bins_DWCE");
	max_dim_x_DUT = iConfig.getParameter<double>("max_dim_x_DUT");
	max_dim_y_DUT = iConfig.getParameter<double>("max_dim_y_DUT");


	TH2F* htmp2;
	TH1F* htmp1;
	std::ostringstream os( std::ostringstream::ate );
	for(int ib = 0; ib<m_NHexaBoards; ib++) {
		os.str("");
		os << "DUT_Occupancy_PerBoard"<<ib;
		htmp2=fs->make<TH2F>(os.str().c_str(), os.str().c_str(), n_bins_DWCE, -max_dim_x_DUT, max_dim_x_DUT, n_bins_DWCE, -max_dim_y_DUT, max_dim_y_DUT);
		m_h_DUT_occupancy_forBoard.insert( std::pair<int,TH2F*>(ib, htmp2) );

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
				os << "DUTOccupancy_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp2 = channel_dir.make<TH2F>("DUTOccupancy", os.str().c_str(), n_bins_DWCE, -max_dim_x_DUT, max_dim_x_DUT, n_bins_DWCE, -max_dim_y_DUT, max_dim_y_DUT);
				h_DUT_occupancy.insert(std::pair<int,TH2F*>(key, htmp2));

				os.str("");
				os << "EnergyVsDUT_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp2=channel_dir.make<TH2F>("RechitEnergyVsDUT", os.str().c_str(), n_bins_DWCE, -max_dim_x_DUT, max_dim_x_DUT, n_bins_DWCE, -max_dim_y_DUT, max_dim_y_DUT);
				m_h_rechitEnergyPerDUT.insert( std::pair<int,TH2F*>(key, htmp2) );
		
				os.str("");
				os << "Energy_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp1=channel_dir.make<TH1F>("RechitEnergy",os.str().c_str(), 33, 1., 100.);
				m_h_rechitEnergy.insert( std::pair<int,TH1F*>(key, htmp1) );


				os.str("");
				os << "EnergyVsDUT_selected_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp2=channel_dir.make<TH2F>("RechitEnergyVsDUT_selected", os.str().c_str(), n_bins_DWCE, -max_dim_x_DUT, max_dim_x_DUT, n_bins_DWCE, -max_dim_y_DUT, max_dim_y_DUT);
				m_h_rechitEnergyPerDUT_selected.insert( std::pair<int,TH2F*>(key, htmp2) );
		
				os.str("");
				os << "Energy_selected_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp1=channel_dir.make<TH1F>("RechitEnergy_selected",os.str().c_str(), 33, 1., 100.);
				m_h_rechitEnergy_selected.insert( std::pair<int,TH1F*>(key, htmp1) );		


				os.str("");
				os << "HGRechitEfficiencyPerDUT"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp2=channel_dir.make<TH2F>("HGRechitEfficiency", os.str().c_str(), n_bins_DWCE, -max_dim_x_DUT, max_dim_x_DUT, n_bins_DWCE, -max_dim_y_DUT, max_dim_y_DUT);
				m_h_rechitEfficiencyPerDUT.insert( std::pair<int,TH2F*>(key, htmp2) );	

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
	int pdgID = rd->pdgID;
	double energy = rd->energy;
	
	if (rd->booleanUserRecords.has("hasDanger")&&rd->booleanUserRecords.get("hasDanger")) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because somthing went wrong"<<std::endl;
		return;
	}

	#ifndef DEBUG
		if (pdgID != 13) {
			std::cout<<"Run is not a dedicated muon run."<<std::endl;
			return;
		}
	#endif

	#ifdef DEBUG
		int eventCounter = rd->event;
		std::cout<<"run: "<<run<<"  energy: "<<energy<<"  pdgID:" << pdgID<<"   eventCounter: "<<eventCounter<<std::endl;
		std::cout<<rd->doubleUserRecords.has("triggerDeltaT_to_TDC")<<"   "<<rd->booleanUserRecords.has("hasValidDWCMeasurement")<<"   "<<rd->booleanUserRecords.has("hasDanger")<<std::endl;
	#endif

	ReadCurrentDWCWindows(run);


	//obtain the track information
	edm::Handle<HGCalTBDWCTrack> dwctrack;
	event.getByToken(DWCTrackToken, dwctrack);
	//Obtain the wire chamber information
	edm::Handle<WireChambers> dwcs;
	event.getByToken(DWCToken, dwcs);
	
	bool vetoEvent = true;
	if (dwctrack->valid) {
		if ((dwctrack->referenceType==15) && (dwctrack->chi2_x<=10.) && (dwctrack->chi2_y<=10.)) { 
			vetoEvent = vetoEvent&&false;
		} else if ((dwctrack->referenceType==13) && (dwctrack->chi2_x<=5.) && (dwctrack->chi2_y<=5.)) { 
			vetoEvent = vetoEvent&&false;
		} else if ((dwctrack->referenceType==14) && (dwctrack->chi2_x<=5.) && (dwctrack->chi2_y<=5.)) { 
			vetoEvent = vetoEvent&&false;
		}
	} else if (dwcs->at(0).goodMeasurement) {
		vetoEvent = vetoEvent&&false;
	} else if (dwcs->at(1).goodMeasurement) {
		vetoEvent = vetoEvent&&false;
	}

	if (vetoEvent) return;

	std::map<int, double> layer_ref_x;
	std::map<int, double> layer_ref_y;

	for(int ib = 0; ib<m_NHexaBoards; ib++) {
		
		int layer = ib+1;		//to do: generalize
  		if (dwctrack->valid) {
	  		if (dwctrack->valid&&(dwctrack->referenceType==15) && (dwctrack->chi2_x<=10.) && (dwctrack->chi2_y<=10.)) { 
	  			layer_ref_x[layer] = dwctrack->DWCExtrapolation_XY(layer).first;
	  			layer_ref_y[layer] = dwctrack->DWCExtrapolation_XY(layer).second;
	  		} else if (dwctrack->valid&&(dwctrack->referenceType==13) && (dwctrack->chi2_x<=5.) && (dwctrack->chi2_y<=5.)) { 
	  			layer_ref_x[layer] = dwctrack->DWCExtrapolation_XY(layer).first;
	  			layer_ref_y[layer] = dwctrack->DWCExtrapolation_XY(layer).second;
	  		} else if (dwctrack->valid&&(dwctrack->referenceType==14) && (dwctrack->chi2_x<=5.) && (dwctrack->chi2_y<=5.)) { 
	  			layer_ref_x[layer] = dwctrack->DWCExtrapolation_XY(layer).first;
	  			layer_ref_y[layer] = dwctrack->DWCExtrapolation_XY(layer).second;
	  		} 
  		} else if (dwcs->at(0).goodMeasurement) {
			layer_ref_x[layer] = dwcs->at(0).x;
			layer_ref_y[layer] = dwcs->at(0).y;
		}
		else if (dwcs->at(1).goodMeasurement) {
			layer_ref_x[layer] = dwcs->at(1).x;
			layer_ref_y[layer] = dwcs->at(1).y;		
		}

		for( size_t _iskiroc=0; _iskiroc<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; _iskiroc++ ) for( size_t _ch=0; _ch<=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; _ch++ ){
			if (_ch%2==1) continue;

			int _key = ib*1000+_iskiroc*100+_ch;

			if ((currentDWCWindows.find(_key) != currentDWCWindows.end()) && (layer_ref_x[layer] > currentDWCWindows[_key][0]) && (layer_ref_x[layer] < currentDWCWindows[_key][1]) && (layer_ref_y[layer] > currentDWCWindows[_key][2]) && (layer_ref_y[layer] < currentDWCWindows[_key][3])) {
				m_h_DUT_occupancy_forBoard[ib]->Fill(layer_ref_x[layer], layer_ref_y[layer], 1.);
				break;
			}
		}
	}

	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);

	
	//fill the rechits:
	for(auto Rechit : *Rechits) {	
		int board = Rechit.board();
		int layer = (Rechit.id()).layer();

		int skiroc = Rechit.skiroc() % 4;
		int channel = Rechit.channel();
  		int key = board*1000+skiroc*100+channel;
  		double energyHG = Rechit.energyHigh();
		
  		double DUT_x = layer_ref_x[layer];
  		double DUT_y = layer_ref_y[layer];


		#ifdef DEBUG
			std::cout<<"Reference type: "<<dwctrack->referenceType<<"   chi2_x: "<<dwctrack->chi2_x<<"   chi2_y: "<<dwctrack->chi2_y<<"    layer: "<<layer<<"  DUT_x: "<<DUT_x<<"   DUT_y: "<<DUT_y<<std::endl;  
		#endif

		h_DUT_occupancy[key]->Fill(DUT_x, DUT_y);



  		m_h_rechitEnergyPerDUT[key]->Fill(DUT_x, DUT_y, energyHG);
  		m_h_rechitEnergy[key]->Fill(energyHG);

		//asses some efficiency for that cell
		//loops over all windows in that board to also assess fake rates
		for( size_t _iskiroc=0; _iskiroc<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; _iskiroc++ ) for( size_t _ch=0; _ch<=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; _ch++ ){
			if (_ch%2==1) continue;
			int _key = board*1000+_iskiroc*100+_ch;
			
			if ((currentDWCWindows.find(_key) != currentDWCWindows.end()) && (DUT_x > currentDWCWindows[_key][0]) && (DUT_x < currentDWCWindows[_key][1]) && (DUT_y > currentDWCWindows[_key][2]) && (DUT_y < currentDWCWindows[_key][3])) {
  				int hit = (energyHG > 0) ? 1: 0;
				m_h_rechitEfficiencyPerDUT[key]->Fill(DUT_x, DUT_y, hit);
				break;
			}
		}


 		if (currentDWCWindows.find(key) == currentDWCWindows.end())	continue;	
		if (DUT_x < currentDWCWindows[key][0] || DUT_x > currentDWCWindows[key][1]) continue;
		if (DUT_y < currentDWCWindows[key][2] || DUT_y > currentDWCWindows[key][3]) continue;

  		m_h_rechitEnergyPerDUT_selected[key]->Fill(DUT_x, DUT_y, energyHG);
  		m_h_rechitEnergy_selected[key]->Fill(energyHG);
	}
	
}// analyze ends here

void MIPFinder::beginJob() {	
	ReadDWCWindows(0);
}

void MIPFinder::endJob() {
	
}

void MIPFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


void MIPFinder::ReadDWCWindows(int fileIndex) {
  	if (fileIndex==(int)pathsToMIPWindowFiles.size()) return;

	std::fstream file; 
	char fragment[100];
	int readCounter = -2;

	WindowMap _parameters;
	  	

	if (pathsToMIPWindowFiles[fileIndex]!=""){
		std::cout<<"Opening: "<<pathsToMIPWindowFiles[fileIndex]<<std::endl;
		file.open(pathsToMIPWindowFiles[fileIndex].c_str(), std::fstream::in);
	}

	int minRun, maxRun;
	if (file.is_open()) {
		file >> fragment;
		minRun = atoi(fragment);
		file >> fragment;
		maxRun = atoi(fragment);
	}

	int iboard = -1, iskiroc = -1, ichannel = -1;
	double DWC_x_min, DWC_x_max, DWC_y_min, DWC_y_max;

	while (file.is_open() && !file.eof()) {		
		if (readCounter!=-2) readCounter++;
			file >> fragment;

		if (std::string(fragment)=="y_max" ) readCounter = -1;  //first parameter is read out

		if (readCounter==0) iboard = atoi(fragment);
		if (readCounter==1) iskiroc = atoi(fragment); 
		if (readCounter==2) ichannel = atoi(fragment);
		if (readCounter==3) DWC_x_min = atof(fragment);
		if (readCounter==4) DWC_x_max = atof(fragment);
		if (readCounter==5) DWC_y_min = atof(fragment);
		if (readCounter==6) { DWC_y_max = atof(fragment);
			int key = iboard*1000+iskiroc*100+ichannel;
			_parameters[key].push_back(DWC_x_min);
			_parameters[key].push_back(DWC_x_max);
			_parameters[key].push_back(DWC_y_min);
			_parameters[key].push_back(DWC_y_max);
			readCounter=-1;
		}
	}
	
	WindowMap ::iterator it;
	for (it=_parameters.begin(); it!=_parameters.end(); it++) {
		std::cout<<"key: "<<it->first;
		for (int i=0; i<4; i++) {
			std::cout<<"  "<<it->second[i];
		}
		std::cout<<std::endl;
	}

	loadedDWCWindows[std::make_pair(minRun, maxRun)] = _parameters;

	return ReadDWCWindows(fileIndex+1);
}


void MIPFinder::ReadCurrentDWCWindows(int this_run) {
	std::map<std::pair<int, int> ,WindowMap  >::iterator it;
	for (it=loadedDWCWindows.begin(); it!=loadedDWCWindows.end(); it++) {
		int run_min = it->first.first;
		int run_max = it->first.second;

		if (this_run>=run_min && (this_run<=run_max || run_max==-1) ) {
			currentDWCWindows = it->second;
			break;
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(MIPFinder);