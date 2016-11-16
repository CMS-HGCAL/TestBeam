/* Write  
 * a description here!
 */

/**
	@Author: Thorben Quast <tquast>
		16 Nov 2016
		thorben.quast@cern.ch / thorben.quast@rwth-aachen.de
*/



// system include files
//#include <memory>
#include <iostream>
#include <string>
//#include "TH2Poly.h"
//#include "TProfile.h"
//#include "TH1F.h"
//#include "TF1.h"
//#include <sstream>
//#include <fstream>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
//#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
//#include "HGCal/Geometry/interface/HGCalTBTopology.h"
//#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
//#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
//#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
//#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
//#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
//#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

//this should go into a config file
double Layer_Z_Positions[16]  = {1.2, 2., 3.5, 4.3, 5.8, 6.3, 8.7, 9.5, 11.4, 12.2, 13.8, 14.6, 16.6, 17.4, 20., 20.8};
int SensorSize = 128;
                    
class Position_Resolution_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

	public:
		explicit Position_Resolution_Analyzer(const edm::ParameterSet&);
		~Position_Resolution_Analyzer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		void analyze(const edm::Event& , const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::EDGetToken HGCalTBRecHitCollection_;
		
		//this is needed to calculate the real x-y positions from u-v
		HGCalTBCellVertices TheCell;
		std::pair<double, double> CellCenterXY;
};



Position_Resolution_Analyzer::Position_Resolution_Analyzer(const edm::ParameterSet& iConfig) {
	// initialization
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
	std::string test_string = iConfig.getParameter<std::string>("randomString");
	std::cout<<"Position_Resolution_Analyzer is initialized with the test_string parameter: "<<test_string<<std::endl;
	
}//constructor ends here

Position_Resolution_Analyzer::~Position_Resolution_Analyzer() {
	return;
}

// ------------ method called for each event  ------------
void Position_Resolution_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {

	//int event_nr = (event.id()).event();
	//double time = event.time().value();

	//opening Rechits
	edm::Handle<HGCalTBRecHitCollection> Rechits;
	event.getByToken(HGCalTBRecHitCollection_, Rechits);

	for(auto Rechit : *Rechits) {
		int layer = (Rechit.id()).layer();
		double iu = (Rechit.id()).iu();
		double iv = (Rechit.id()).iv();
		//double energyLow = Rechit.energyLow();
		//double energyHigh = Rechit.energyHigh();
		double energy = Rechit.energy();

		CellCenterXY = TheCell.GetCellCentreCoordinatesForPlots((Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv(), SensorSize);
		double iux = CellCenterXY.first;
		double ivy = CellCenterXY.second;
		std::cout<<layer<<"  "<<iu<<"  "<<iv<<"  "<<iux<<"  "<<ivy<<"  "<<energy<<std::endl;
		

		//(Rechit.id()).layer(), (Rechit.id()).sensorIU(), (Rechit.id()).sensorIV(), (Rechit.id()).iu(), (Rechit.id()).iv()
		//Rechit.energyHigh(), Rechit.energyLow()
	}
	
}// analyze ends here



void Position_Resolution_Analyzer::beginJob() {
}


void Position_Resolution_Analyzer::endJob() {
}


void
Position_Resolution_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Position_Resolution_Analyzer);
