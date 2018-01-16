/* 
 * Writing of ROOT Trees given the user records.
 */

/**
	@Author: Thorben Quast <tquast>
		22 November 2017
		thorben.quast@cern.ch / thorben.quast@rwth-aachen.de
*/


// system include files
#include <iostream>
#include <string>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HGCal/DataFormats/interface/HGCalTBRunData.h"	//for the runData type definition
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "TFile.h"
#include "TTree.h"

//#define DEBUG

class NTupelizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit NTupelizer(const edm::ParameterSet&);
		~NTupelizer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		void analyze(const edm::Event& , const edm::EventSetup&) override;

		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		edm::EDGetTokenT<UserRecords<double>> UserRecordToken;
		std::vector<std::string> UserRecordKeys;
		size_t N_keys;

		TTree* outTree;

		std::vector<double> dataToFill;
};

NTupelizer::NTupelizer(const edm::ParameterSet& iConfig) {	
	UserRecordToken = consumes<UserRecords<double> >(iConfig.getParameter<edm::InputTag>("USERRECORDS"));			
	UserRecordKeys = iConfig.getParameter<std::vector<std::string> >("UserRecordKeys");
	N_keys = UserRecordKeys.size();

	usesResource("TFileService");
	outTree = fs->make<TTree>("variables", "variables");
	for (size_t i=0; i<N_keys; i++) dataToFill.push_back(-1.);

	for (size_t i=0; i<N_keys; i++) 
		outTree->Branch(UserRecordKeys[i].c_str(), &dataToFill[i], (UserRecordKeys[i]+"/D").c_str());
}

NTupelizer::~NTupelizer() {
	return;
}

// ------------ method called for each event  ------------
void NTupelizer::analyze(const edm::Event& event, const edm::EventSetup& setup) {

	edm::Handle<UserRecords<double> > ur;
	event.getByToken(UserRecordToken, ur);
	
	for (size_t i=0; i<N_keys; i++) {
		dataToFill[i] = -999.;
		if (!ur->has(UserRecordKeys[i])) continue;
		dataToFill[i]=1.0*ur->get(UserRecordKeys[i]);
	}
	outTree->Fill();

	#ifdef DEBUG
		std::cout<<ur->has("eventID")<<"  "<<ur->get("eventID")<<std::endl;
		std::cout<<ur->has("run")<<"  "<<ur->get("run")<<std::endl;
		std::cout<<ur->has("pdgID")<<"  "<<ur->get("pdgID")<<std::endl;
		std::cout<<ur->has("beamEnergy")<<"  "<<ur->get("beamEnergy")<<std::endl;
	#endif
	
}// analyze ends here



void NTupelizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NTupelizer);