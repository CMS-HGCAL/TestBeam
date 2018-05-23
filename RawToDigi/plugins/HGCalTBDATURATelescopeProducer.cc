#include "HGCal/RawToDigi/plugins/HGCalTBDATURATelescopeProducer.h"

//#define DEBUG

HGCalTBDATURATelescopeProducer::HGCalTBDATURATelescopeProducer(const edm::ParameterSet& cfg) {
	RunDataToken = consumes<RunData>(cfg.getParameter<edm::InputTag>("RUNDATA"));
	inputFile = cfg.getParameter<std::string>("inputFile");
    outputCollectionName = cfg.getParameter<std::string>("OutputCollectionName");
    SkipFirstNEvents = cfg.getParameter<int>("SkipFirstNEventsInTelescopeFile");
    m_layerPositionFile = cfg.getParameter<std::string>("layerPositionFile");
    m_PIStagePositionFile = cfg.getParameter<std::string>("PIStagePositionFile");


    produces<std::vector<HGCalTBDATURATelescopeData> >(outputCollectionName);
    produces<RunData>("FullRunData");
}

void HGCalTBDATURATelescopeProducer::beginJob() {
    rootFile = new TFile(inputFile.c_str(), "READ");
    tree = (TTree*)rootFile->Get("corryvreckan/HGCalTBDataOutput/trackClusters");           //hard coded value, might be subject to change
    
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


    //reads the layer positions

    std::fstream file; 
    char fragment[100];
    int readCounter = -1;

    file.open(m_layerPositionFile.c_str(), std::fstream::in);
    std::cout<<"Reading file "<<m_layerPositionFile<<" -open: "<<file.is_open()<<std::endl;
    int layer=0;
    while (file.is_open() && !file.eof()) {
        readCounter++;
        file >> fragment;
        if (readCounter==0) layer=atoi(fragment);
        if (readCounter==1) {
            layerPositions[layer]=atof(fragment)/10.;       //values are given in mm and should be converted into cm
            readCounter=-1;
        }
    }
    file.close();

    //reads the PI stage positions
    readCounter = -1;
    file.open(m_PIStagePositionFile.c_str(), std::fstream::in);
    std::cout<<"Reading file "<<m_PIStagePositionFile<<" -open: "<<file.is_open()<<std::endl;
    int run=0;
    float pos_x = 0;
    float pos_y = 0;
    while (file.is_open() && !file.eof()) {
        readCounter++;
        file >> fragment;
        if (readCounter==0) continue;
        if (readCounter==1) run=atoi(fragment);
        if (readCounter==2) continue;
        if (readCounter==3) pos_x=atof(fragment);
        if (readCounter==4) continue;
        if (readCounter==5) {
            pos_y=atof(fragment);
            PIStagePositions[run]=std::make_pair(pos_x, pos_y);
            readCounter=-1;
        }
    }

}

void HGCalTBDATURATelescopeProducer::produce(edm::Event& event, const edm::EventSetup& iSetup) {

 	//get the relevant event information
	edm::Handle<RunData> rd;
	event.getByToken(RunDataToken, rd);

    std::unique_ptr<std::vector<HGCalTBDATURATelescopeData> > DATURATracks(new std::vector<HGCalTBDATURATelescopeData>);

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
    
    if (PIStagePositions.find((int)rd->run) != PIStagePositions.end()) {
        rd_full->doubleUserRecords.add("PIStagePosition_X", PIStagePositions[(int)rd->run].first);
        rd_full->doubleUserRecords.add("PIStagePosition_Y", PIStagePositions[(int)rd->run].second);
    }
        
    
    if (tree_Ntracks<1) {
        rd_full->booleanUserRecords.add("hasValidDATURAMeasurement", false);   
    }
    else {
        for (int ntrack=1; ntrack <= tree_Ntracks; ntrack++) {
            
           HGCalTBDATURATelescopeData DATURATelescopeTrack(ntrack);
           LineFitter TripletTrack1X;
           LineFitter TripletTrack1Y;
           LineFitter TripletTrack2X;
           LineFitter TripletTrack2Y;           
            for (int MIMOSA_index=1; MIMOSA_index<=3; MIMOSA_index++) {
                DATURATelescopeTrack.addPointForTracking(tree_clusterX[MIMOSA_index]->at(ntrack-1), tree_clusterY[MIMOSA_index]->at(ntrack-1), tree_clusterZ[MIMOSA_index]->at(ntrack-1), 0.0184, 0.0184);
                TripletTrack1X.addPoint(tree_clusterZ[MIMOSA_index]->at(ntrack-1), tree_clusterX[MIMOSA_index]->at(ntrack-1), 0.0184);
                TripletTrack1Y.addPoint(tree_clusterZ[MIMOSA_index]->at(ntrack-1), tree_clusterY[MIMOSA_index]->at(ntrack-1), 0.0184);
            }   
            for (int MIMOSA_index=4; MIMOSA_index<=6; MIMOSA_index++) {
                DATURATelescopeTrack.addPointForTracking(tree_clusterX[MIMOSA_index]->at(ntrack-1), tree_clusterY[MIMOSA_index]->at(ntrack-1), tree_clusterZ[MIMOSA_index]->at(ntrack-1), 0.0184, 0.0184);
                TripletTrack2X.addPoint(tree_clusterZ[MIMOSA_index]->at(ntrack-1), tree_clusterX[MIMOSA_index]->at(ntrack-1), 0.0184);
                TripletTrack2Y.addPoint(tree_clusterZ[MIMOSA_index]->at(ntrack-1), tree_clusterY[MIMOSA_index]->at(ntrack-1), 0.0184);
            }   

            TripletTrack1X.fit();       TripletTrack1Y.fit();
            TripletTrack2X.fit();       TripletTrack2Y.fit();

            for (std::map<int, double>::iterator layerIt=layerPositions.begin(); layerIt!=layerPositions.end(); layerIt++) {
                double layer_ref_x = TripletTrack2X.eval(layerIt->second);
                double layer_ref_y = TripletTrack2Y.eval(layerIt->second);
                double layer_ref_x_chi2 = TripletTrack2X.GetChisquare();
                double layer_ref_y_chi2 = TripletTrack2Y.GetChisquare();
                
                //mean between two triplets for the intermediate layer
                if (layerIt->first <=1 ) {
                    layer_ref_x += TripletTrack1X.eval(layerIt->second);
                    layer_ref_y += TripletTrack1Y.eval(layerIt->second);
                    layer_ref_x_chi2 += TripletTrack1X.GetChisquare();
                    layer_ref_y_chi2 += TripletTrack1Y.GetChisquare();

                    layer_ref_x /= 2.;
                    layer_ref_y /= 2.;
                    layer_ref_x_chi2 /= 2.;
                    layer_ref_y_chi2 /= 2.;
                }
                DATURATelescopeTrack.addLayerReference(layerIt->first, layer_ref_x, layer_ref_y, layer_ref_x_chi2, layer_ref_y_chi2);
            }

            float kinkAngleX_DUT1 = atan(TripletTrack1X.getM())-atan(TripletTrack2X.getM());
            float kinkAngleY_DUT1 = atan(TripletTrack1Y.getM())-atan(TripletTrack2Y.getM());
            DATURATelescopeTrack.floatUserRecords.add("kinkAngleX_DUT1", kinkAngleX_DUT1);
            DATURATelescopeTrack.floatUserRecords.add("kinkAngleY_DUT1", kinkAngleY_DUT1);
            DATURATracks->push_back(DATURATelescopeTrack);
        }

        #ifdef DEBUG
            for (size_t i=0; i<DATURATracks->size(); i++) std::cout<<"Adding a telescope track with "<<DATURATracks->at(i).Extrapolation_XY(1).first<<", "<<DATURATracks->at(i).Extrapolation_XY(1).second<<"   ("<<DATURATracks->at(i).chi2_x<<", "<<DATURATracks->at(i).chi2_y<<")"<<std::endl;
        #endif
        rd_full->booleanUserRecords.add("hasValidDATURAMeasurement", true);
    }
    event.put(std::move(rd_full), "FullRunData");     
    event.put(std::move(DATURATracks), outputCollectionName);

}


void HGCalTBDATURATelescopeProducer::endJob() {
    delete tree;
    delete rootFile;
}

DEFINE_FWK_MODULE(HGCalTBDATURATelescopeProducer);
