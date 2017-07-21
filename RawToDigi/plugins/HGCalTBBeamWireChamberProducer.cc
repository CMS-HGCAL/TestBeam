#include "HGCal/RawToDigi/plugins/HGCalTBBeamWireChamberProducer.h"

//#define DEBUG

HGCalTBBeamWireChamberProducer::HGCalTBBeamWireChamberProducer(const edm::ParameterSet& cfg) {
	RunDataToken = consumes<RunData>(cfg.getParameter<edm::InputTag>("RUNDATA"));
	MWCToken = consumes<WireChambers>(cfg.getParameter<edm::InputTag>("MWCHAMBERS"));		//for cross-check
	inputFile = cfg.getParameter<std::string>("inputFile");
}

void HGCalTBBeamWireChamberProducer::beginJob() {
    rootFile = new TFile(inputFile.c_str(), "READ");
    tree = (TTree*)rootFile->Get("dwc_ntupelizer/dwc_reco");
    
	tree->SetBranchAddress("run", &run, &b_run);
    tree->SetBranchAddress("event", &eventId ,&b_event);   
    tree->SetBranchAddress("goodDWC_Measurement", &goodDWC_Measurement ,&b_goodDWC_Measurement);   
    tree->SetBranchAddress("reco1_x", &reco1_x ,&b_reco1_x);   
    tree->SetBranchAddress("reco1_y", &reco1_y ,&b_reco1_y);   
    tree->SetBranchAddress("z1", &z1 ,&b_z1);   
    tree->SetBranchAddress("reco2_x", &reco2_x ,&b_reco2_x);   
    tree->SetBranchAddress("reco2_y", &reco2_y ,&b_reco2_y);   
    tree->SetBranchAddress("z2", &z2 ,&b_z2);   
    tree->SetBranchAddress("reco3_x", &reco3_x ,&b_reco3_x);   
    tree->SetBranchAddress("reco3_y", &reco3_y ,&b_reco3_y);   
    tree->SetBranchAddress("z3", &z3 ,&b_z3);   
    tree->SetBranchAddress("reco4_x", &reco4_x ,&b_reco4_x);   
    tree->SetBranchAddress("reco4_y", &reco4_y ,&b_reco4_y);   
    tree->SetBranchAddress("z4", &z4 ,&b_z4);     

    loaded_run = -1;
}

void HGCalTBBeamWireChamberProducer::produce(edm::Event& event, const edm::EventSetup& iSetup) {

 	//get the relevant event information
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);
	if (rd->event == -1) return;

	//get the multi wire chambers
	edm::Handle<WireChambers> dwcs;
	event.getByToken(MWCToken, dwcs);

	if (rd->run != loaded_run) {
		loadRun(rd->run);
	}

	//cross-check
	#ifdef DEBUG
		std::cout<<"Event: "<<rd->event<<std::endl;
		std::cout<<"reco1_x: "<<dwcs->at(0).x<<" vs. "<<reco1_x_loaded[rd->event]<<std::endl;
		std::cout<<"reco1_y: "<<dwcs->at(0).y<<" vs. "<<reco1_y_loaded[rd->event]<<std::endl;
		std::cout<<"z1: "<<dwcs->at(0).z<<" vs. "<<z1_loaded[rd->event]<<std::endl;

		std::cout<<"reco2_x: "<<dwcs->at(1).x<<" vs. "<<reco2_x_loaded[rd->event]<<std::endl;
		std::cout<<"reco2_y: "<<dwcs->at(1).y<<" vs. "<<reco2_y_loaded[rd->event]<<std::endl;
		std::cout<<"z2: "<<dwcs->at(1).z<<" vs. "<<z2_loaded[rd->event]<<std::endl;

		std::cout<<"reco3_x: "<<dwcs->at(2).x<<" vs. "<<reco3_x_loaded[rd->event]<<std::endl;
		std::cout<<"reco3_y: "<<dwcs->at(2).y<<" vs. "<<reco3_y_loaded[rd->event]<<std::endl;
		std::cout<<"z3: "<<dwcs->at(2).z<<" vs. "<<z3_loaded[rd->event]<<std::endl;

		std::cout<<"reco4_x: "<<dwcs->at(3).x<<" vs. "<<reco4_x_loaded[rd->event]<<std::endl;
		std::cout<<"reco4_y: "<<dwcs->at(3).y<<" vs. "<<reco4_y_loaded[rd->event]<<std::endl;
		std::cout<<"z4: "<<dwcs->at(3).z<<" vs. "<<z4_loaded[rd->event]<<std::endl<<std::endl;
	#endif

}

void HGCalTBBeamWireChamberProducer::loadRun(int loading_run) {
	#ifdef DEBUG
    	std::cout<<"Clearing run "<<loaded_run<<std::endl;
	#endif
    reco1_x_loaded.clear(); reco1_y_loaded.clear(); z1_loaded.clear();
    reco2_x_loaded.clear(); reco2_y_loaded.clear(); z2_loaded.clear();
    reco3_x_loaded.clear(); reco3_y_loaded.clear(); z3_loaded.clear();
    reco4_x_loaded.clear(); reco4_y_loaded.clear(); z4_loaded.clear();	

	#ifdef DEBUG
    	std::cout<<"Loading run "<<loading_run<<std::endl;
    #endif
    for (size_t i=0; i<(size_t) tree->GetEntries(); i++){
    	tree->GetEntry(i);
    	if (run!=loading_run) continue;
    	reco1_x_loaded[eventId] = reco1_x;
    	reco1_y_loaded[eventId] = reco1_y;
    	z1_loaded[eventId] = z1;
      	reco2_x_loaded[eventId] = reco2_x;
    	reco2_y_loaded[eventId] = reco2_y;
    	z2_loaded[eventId] = z2;
       	reco3_x_loaded[eventId] = reco3_x;
    	reco3_y_loaded[eventId] = reco3_y;
    	z3_loaded[eventId] = z3;
       	reco4_x_loaded[eventId] = reco4_x;
    	reco4_y_loaded[eventId] = reco4_y;
    	z4_loaded[eventId] = z4;
    }
	#ifdef DEBUG
    	std::cout<<"Loaded run "<<loading_run<<std::endl;
    #endif

	loaded_run = loading_run;
}

void HGCalTBBeamWireChamberProducer::endJob() {
    delete tree;
    delete rootFile;

    reco1_x_loaded.clear(); reco1_y_loaded.clear(); z1_loaded.clear();
    reco2_x_loaded.clear(); reco2_y_loaded.clear(); z2_loaded.clear();
    reco3_x_loaded.clear(); reco3_y_loaded.clear(); z3_loaded.clear();
    reco4_x_loaded.clear(); reco4_y_loaded.clear(); z4_loaded.clear();	
}

DEFINE_FWK_MODULE(HGCalTBBeamWireChamberProducer);
