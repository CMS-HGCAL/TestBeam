#include "HGCal/RawToDigi/plugins/HGCalTBDATURATelescopeProducer.h"

//#define DEBUG

HGCalTBDATURATelescopeProducer::HGCalTBDATURATelescopeProducer(const edm::ParameterSet& cfg) {
	RunDataToken = consumes<RunData>(cfg.getParameter<edm::InputTag>("RUNDATA"));
	inputFile = cfg.getParameter<std::string>("inputFile");
    outputCollectionName = cfg.getParameter<std::string>("OutputCollectionName");
    SkipFirstNEvents = cfg.getParameter<int>("SkipFirstNEventsInTelescopeFile");


    produces<std::map<int, WireChamberData> >(outputCollectionName);
    produces<RunData>("FullRunData");
}

void HGCalTBDATURATelescopeProducer::beginJob() {
    rootFile = new TFile(inputFile.c_str(), "READ");
    tree = (TTree*)rootFile->Get("corryvreckan/HGCalTBDataOutput/trackClusters");
    
    tree->SetBranchAddress("EventID", &tree_eventID, &b_event);
    tree->SetBranchAddress("chi2", &tree_track_chi2, &b_chi2);
	tree->SetBranchAddress("NTracks", &tree_Ntracks, &b_Ntracks);

    for (int MIMOSA_index=1; MIMOSA_index<=6; MIMOSA_index++) {
        b_clusterX[MIMOSA_index] = new TBranch;
        b_clusterY[MIMOSA_index] = new TBranch;
        b_clusterZ[MIMOSA_index] = new TBranch;
        b_absorber[MIMOSA_index] = new TBranch;


        tree->SetBranchAddress(("associatedClusterX_MIMOSA26_" + std::to_string(MIMOSA_index)).c_str(), &tree_clusterX[MIMOSA_index], &b_clusterX.at(MIMOSA_index));
        tree->SetBranchAddress(("associatedClusterY_MIMOSA26_" + std::to_string(MIMOSA_index)).c_str(), &tree_clusterY[MIMOSA_index], &b_clusterY.at(MIMOSA_index));
        tree->SetBranchAddress(("associatedClusterZ_MIMOSA26_" + std::to_string(MIMOSA_index)).c_str(), &tree_clusterZ[MIMOSA_index], &b_clusterZ.at(MIMOSA_index));
        tree->SetBranchAddress(("associatedAbsorber_MIMOSA26_" + std::to_string(MIMOSA_index)).c_str(), &tree_absorber[MIMOSA_index], &b_absorber.at(MIMOSA_index));

        tree_clusterX[MIMOSA_index] = 0;
        tree_clusterY[MIMOSA_index] = 0;
        tree_clusterZ[MIMOSA_index] = 0;
        tree_absorber[MIMOSA_index] = 0;

    }

    tree_track_chi2=0;      //set pointer to vector to 0
    
}

void HGCalTBDATURATelescopeProducer::produce(edm::Event& event, const edm::EventSetup& iSetup) {

 	//get the relevant event information
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);

    std::unique_ptr<std::map<int, WireChamberData> > TelescopePlanes(new std::map<int, WireChamberData>);

    tree->GetEntry(rd->event - 1 + SkipFirstNEvents);
    #ifdef DEBUG
        std::cout<<"Run event: "<<rd->event<<"  vs. eventID in tree: "<<tree_eventID<<"   Number of tracks:" << tree_Ntracks<<std::endl;
        for (int nt=0; nt<tree_Ntracks; nt++) {   
                std::cout<<"Track "<<nt+1<<": "<<tree_track_chi2->at(nt)<<std::endl;
            for (int MIMOSA_index=1; MIMOSA_index<=6; MIMOSA_index++) {
                std::cout<<tree_clusterX[MIMOSA_index]->at(nt)<<"  ,  "<<tree_clusterY[MIMOSA_index]->at(nt)<<"  ,  "<<tree_clusterZ[MIMOSA_index]->at(nt)<<std::endl;
            }
        }
    #endif

    std::unique_ptr<RunData> rd_full(new RunData);
    rd_full->configuration = rd->configuration;
    rd_full->run = rd->run;
    rd_full->trigger = rd->trigger;
    rd_full->event = rd->event;
    rd_full->energy = rd->energy;
    rd_full->runType = rd->runType;
    rd_full->pdgID = rd->pdgID;
    if (rd->booleanUserRecords.has("hasDanger")) rd_full->booleanUserRecords.add("hasDanger", rd->booleanUserRecords.get("hasDanger"));
    if (rd->doubleUserRecords.has("trueEnergy")) rd_full->doubleUserRecords.add("trueEnergy", rd->doubleUserRecords.get("trueEnergy"));

    if (tree_Ntracks!=1) {
        rd_full->booleanUserRecords.add("hasValidDWCMeasurement", false);
        
    }
    else {
        rd_full->booleanUserRecords.add("hasValidDWCMeasurement", true);
        //just take the clusters with the best chi2
        for (int MIMOSA_index=1; MIMOSA_index<=6; MIMOSA_index++) {
            WireChamberData* TelescopePlane = new WireChamberData();

            TelescopePlane->ID = MIMOSA_index;
            TelescopePlane->goodMeasurement_X = true;
            TelescopePlane->goodMeasurement_Y = true;
            TelescopePlane->goodMeasurement = (TelescopePlane->goodMeasurement_X && TelescopePlane->goodMeasurement_Y);
            TelescopePlane->x = tree_clusterX[MIMOSA_index]->at(0);
            TelescopePlane->res_x = 0.0184; //18.4microns pitch of the MIMOSA telescope planes
            TelescopePlane->y = tree_clusterY[MIMOSA_index]->at(0);
            TelescopePlane->res_y = 0.0184; //18.4microns pitch of the MIMOSA telescope planes
            TelescopePlane->z = tree_clusterZ[MIMOSA_index]->at(0);
            TelescopePlane->averageHitMultiplicty = 1;                    

            (*TelescopePlanes)[MIMOSA_index-1] = *TelescopePlane;
        }   
    }
    event.put(std::move(rd_full), "FullRunData");     
    event.put(std::move(TelescopePlanes), outputCollectionName);

}


void HGCalTBDATURATelescopeProducer::endJob() {
    delete tree;
    delete rootFile;
}

DEFINE_FWK_MODULE(HGCalTBDATURATelescopeProducer);
