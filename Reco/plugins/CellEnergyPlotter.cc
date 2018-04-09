/* 
 * Cell energy plotting.
 */

/**
	@Author: Thorben Quast <tquast>
		08 April 2018
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
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBCommonModeNoise.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"  
#include <sstream>
#include <fstream>
#include <iomanip>
#include <set>


bool rejectFromCommonModeNoise(edm::Handle<std::map<int, commonModeNoise> > &cmMap, int iski, int criterion=0) {
	if (criterion==0) return false;
	else if (criterion==1) {
		for (size_t ts=1; ts<=6; ts++)	if (fabs(cmMap->at(iski).fullHG[ts]) > 200.) return true;
		return false;
	}
	else if (criterion==2) {
		for (size_t ts=1; ts<=6; ts++)	if (cmMap->at(iski).fullHG[ts] < -200.) return true;
		return false;
	}
	else return false;

}

class CellEnergyPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit CellEnergyPlotter(const edm::ParameterSet&);
		~CellEnergyPlotter();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<HGCalTBRecHitCollection> HGCalTBRecHitCollection_Token;	 	
		edm::EDGetTokenT<std::map<int, commonModeNoise> > CommonModeNoiseMap_Token;	 	
		edm::EDGetTokenT<RunData> RunDataToken;	
	
		std::string m_electronicMap;
		std::string m_detectorLayoutFile;
		struct {
			HGCalElectronicsMap emap_;
			HGCalTBDetectorLayout layout_;
		} essource_;

		int m_NHexaBoards;

		std::map<int,TH1F*> m_h_rechitEnergyHG;
		std::map<int,TH1F*> m_h_rechitEnergyLG;
		std::map<int,TH1F*> m_h_rechitEnergy;
		std::map<int,TH1F*> m_h_Occupancy;		//based on successful rechit energy reco
		std::map<int,TH2F*> m_h_rechitEnergyHGvsChannel;
		std::map<int,TH2F*> m_h_rechitEnergyLGvsChannel;
		std::map<int,TH2F*> m_h_rechitEnergyvsChannel;

		int commonModeNoiseRejectionType;
		bool rejectFromCommonModeNoise(edm::Handle<std::map<int, commonModeNoise> > &cmMap, int iski, int criterion);
};

CellEnergyPlotter::CellEnergyPlotter(const edm::ParameterSet& iConfig) {	
	
	usesResource("TFileService");
	HGCalTBRecHitCollection_Token = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	CommonModeNoiseMap_Token = consumes<std::map<int, commonModeNoise>>(iConfig.getParameter<edm::InputTag>("HGCALTBCOMMONMODENOISE"));
	RunDataToken= consumes<RunData>(iConfig.getParameter<edm::InputTag>("RUNDATA"));
	
	m_electronicMap = iConfig.getUntrackedParameter<std::string>("ElectronicMap","HGCal/CondObjects/data/map_CERN_Hexaboard_28Layers_AllFlipped.txt");
	m_detectorLayoutFile = iConfig.getUntrackedParameter<std::string>("DetectorLayout","HGCal/CondObjects/data/layerGeom_oct2017_h2_17layers.txt");
	m_NHexaBoards= iConfig.getUntrackedParameter<int>("NHexaBoards", 10);


	commonModeNoiseRejectionType = iConfig.getParameter<int>("commonModeNoiseRejectionType");


	TH1F* htmp1;
	TH2F* htmp2;
	std::ostringstream os( std::ostringstream::ate );
	for(int ib = 0; ib<m_NHexaBoards; ib++) {

		for( size_t iski=0; iski<HGCAL_TB_GEOMETRY::N_SKIROC_PER_HEXA; iski++ ){
			os.str("");
			os<<"HexaBoard"<<ib<<"_Skiroc"<<iski;
			TFileDirectory chip_dir = fs->mkdir( os.str().c_str() );
			int keyChip = ib*1000+iski*100;

			os.str("");
			os << "HG_board_"<<ib<<"_chip_"<<iski;
			htmp2=chip_dir.make<TH2F>("RechitEnergyHGVsChannel",os.str().c_str(), 64., -0.5, 63.5, 50, 1., 150.);
			m_h_rechitEnergyHGvsChannel.insert( std::pair<int,TH2F*>(keyChip, htmp2) );

			os.str("");
			os << "LG_board_"<<ib<<"_chip_"<<iski;
			htmp2=chip_dir.make<TH2F>("RechitEnergyLGVsChannel",os.str().c_str(), 64., -0.5, 63.5, 50, 1., 20.);
			m_h_rechitEnergyLGvsChannel.insert( std::pair<int,TH2F*>(keyChip, htmp2) );

			os.str("");
			os << "Energy_board_"<<ib<<"_chip_"<<iski;
			htmp2=chip_dir.make<TH2F>("RechitEnergyVsChannel",os.str().c_str(), 64., -0.5, 63.5, 50, 0., 4.);
			m_h_rechitEnergyvsChannel.insert( std::pair<int,TH2F*>(keyChip, htmp2) );

			os.str("");
			os << "Occupancy_board"<<ib<<"_chip_"<<iski;
			htmp1=chip_dir.make<TH1F>("Occupancy",os.str().c_str(), 64., -0.5, 63.5);
			m_h_Occupancy.insert( std::pair<int,TH1F*>(keyChip, htmp1) );

			for( size_t ichan=0; ichan<=HGCAL_TB_GEOMETRY::N_CHANNELS_PER_SKIROC; ichan++ ){
				if (ichan%2==1) continue;
				int keyChannel = ib*1000+iski*100+ichan;
				os.str("");
				os << "Channel" << ichan;
				TFileDirectory channel_dir = chip_dir.mkdir( os.str().c_str() );

				os.str("");
				os << "HG_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp1=channel_dir.make<TH1F>("RechitEnergyHG",os.str().c_str(), 50, 1., 150.);
				m_h_rechitEnergyHG.insert( std::pair<int,TH1F*>(keyChannel, htmp1) );

				os.str("");
				os << "LG_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp1=channel_dir.make<TH1F>("RechitEnergyLG",os.str().c_str(), 50, 1., 20.);
				m_h_rechitEnergyLG.insert( std::pair<int,TH1F*>(keyChannel, htmp1) );

				os.str("");
				os << "Energy_board_"<<ib<<"_chip_"<<iski<<"_channel_"<<ichan;
				htmp1=channel_dir.make<TH1F>("RechitEnergy",os.str().c_str(), 50, 0., 4.);
				m_h_rechitEnergy.insert( std::pair<int,TH1F*>(keyChannel, htmp1) );

			}			
		}
	}


}

CellEnergyPlotter::~CellEnergyPlotter() {
	return;
}

// ------------ method called for each event  ------------
void CellEnergyPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup) {

	edm::Handle<RunData> rd;
 	//get the relevant event information
	event.getByToken(RunDataToken, rd);
	
	int evId = event.id().event();
	int run = rd->run;
	//int pdgID = rd->pdgID;
	double energy = rd->energy;
	
	if (rd->booleanUserRecords.has("hasDanger")&&rd->booleanUserRecords.get("hasDanger")) {
		std::cout<<"Event "<<evId<<" of run "<<run<<" ("<<energy<<"GeV)  is skipped because somthing went wrong"<<std::endl;
		return;
	}


	edm::Handle<std::map<int, commonModeNoise>> cmMap;
	if (commonModeNoiseRejectionType) {
		event.getByToken(CommonModeNoiseMap_Token, cmMap);
	}

	
	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_Token, Rechits);

	
	//fill the rechits:
	for(auto Rechit : *Rechits) {	
		//int layer = (Rechit.id()).layer();

		HGCalTBElectronicsId eid( essource_.emap_.detId2eid( Rechit.id().rawId() ) );
		int skiroc = eid.iskiroc_rawhit();
		int board = skiroc / 4;
		if (commonModeNoiseRejectionType&&rejectFromCommonModeNoise(cmMap, skiroc, commonModeNoiseRejectionType)) continue;
		
		skiroc = eid.iskiroc_rawhit() % 4;
		int channel = eid.ichan();

		int keyChip = board*1000+skiroc*100;
  		int keyChannel = board*1000+skiroc*100+channel;
  		
  		if (!Rechit.checkFlag(HGCalTBRecHit::kGood)) continue;
  		double energyHG = (!Rechit.checkFlag(HGCalTBRecHit::kHighGainSaturated)) ? Rechit.energyHigh() : -1;
  		double energyLG = (!Rechit.checkFlag(HGCalTBRecHit::kLowGainSaturated)) ? Rechit.energyLow() : -1;
  		double energy = Rechit.energy();
		
		m_h_Occupancy[keyChip]->Fill(channel);		
		m_h_rechitEnergyHGvsChannel[keyChip]->Fill(channel, energyHG);
		m_h_rechitEnergyLGvsChannel[keyChip]->Fill(channel, energyLG);
		m_h_rechitEnergyvsChannel[keyChip]->Fill(channel, energy);
		m_h_rechitEnergyHG[keyChannel]->Fill(energyHG);
		m_h_rechitEnergyLG[keyChannel]->Fill(energyLG);
		m_h_rechitEnergy[keyChannel]->Fill(energy);
	}
	
}// analyze ends here

void CellEnergyPlotter::beginJob() {	

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

void CellEnergyPlotter::endJob() {
	
}

void CellEnergyPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


bool CellEnergyPlotter::rejectFromCommonModeNoise(edm::Handle<std::map<int, commonModeNoise> > &cmMap, int iski, int criterion=0) {
	if (criterion==0) return false;
	else if (criterion==1) {
		for (size_t ts=1; ts<=6; ts++)	if (fabs(cmMap->at(iski).fullHG[ts]) > 200.) return true;
		return false;
	}
	else if (criterion==2) {
		for (size_t ts=1; ts<=6; ts++)	if (cmMap->at(iski).fullHG[ts] < -200.) return true;
		return false;
	}
	else return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CellEnergyPlotter);