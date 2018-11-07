#include "HGCal/RawToDigi/plugins/HGCalTBCAENProducer.h"
#include <fstream>

//#define DEBUG

HGCalTBBeamWireChamberProducer::HGCalTBBeamWireChamberProducer(const edm::ParameterSet& cfg) {
    RunDataToken = consumes<RunData>(cfg.getParameter<edm::InputTag>("RUNDATA"));
    inputFile = cfg.getParameter<std::string>("inputFile");
    outputCollectionName = cfg.getParameter<std::string>("OutputCollectionName");

    produces<std::map<int, WireChamberData> >(outputCollectionName);
    produces<RunData>("FullRunData");

    tree = NULL;
}

void HGCalTBBeamWireChamberProducer::beginJob() {
    std::ifstream infile(inputFile);
    if (!infile.good()) {
        rootFile = NULL; 
        tree = NULL;
        return;
    }

    rootFile = new TFile(inputFile.c_str(), "READ");
    tree = (TTree*)rootFile->Get("dwc_ntupelizer/dwc_reco");

    tree->SetBranchAddress("run", &run, &b_run);
    tree->SetBranchAddress("event", &eventId , &b_event);
    tree->SetBranchAddress("goodDWC_Measurement", &goodDWC_Measurement , &b_goodDWC_Measurement);
    tree->SetBranchAddress("triggerTimeDifference", &triggerTimeDiff , &b_triggerTimeDiff);
    tree->SetBranchAddress("reco1_x", &reco1_x , &b_reco1_x);
    tree->SetBranchAddress("res1_x", &res1_x , &b_res1_x);
    tree->SetBranchAddress("reco1_y", &reco1_y , &b_reco1_y);
    tree->SetBranchAddress("res1_y", &res1_y , &b_res1_y);
    tree->SetBranchAddress("z1", &z1 , &b_z1);
    tree->SetBranchAddress("dwc1_multiplicity", &averageHitMultiplicty1 , &b_averageHitMultiplicty1);
    tree->SetBranchAddress("reco2_x", &reco2_x , &b_reco2_x);
    tree->SetBranchAddress("res2_x", &res2_x , &b_res2_x);
    tree->SetBranchAddress("reco2_y", &reco2_y , &b_reco2_y);
    tree->SetBranchAddress("res2_y", &res2_y , &b_res2_y);
    tree->SetBranchAddress("z2", &z2 , &b_z2);
    tree->SetBranchAddress("dwc2_multiplicity", &averageHitMultiplicty2 , &b_averageHitMultiplicty2);
    tree->SetBranchAddress("reco3_x", &reco3_x , &b_reco3_x);
    tree->SetBranchAddress("res3_x", &res3_x , &b_res3_x);
    tree->SetBranchAddress("reco3_y", &reco3_y , &b_reco3_y);
    tree->SetBranchAddress("res3_y", &res3_y , &b_res3_y);
    tree->SetBranchAddress("z3", &z3 , &b_z3);
    tree->SetBranchAddress("dwc3_multiplicity", &averageHitMultiplicty3 , &b_averageHitMultiplicty3);
    tree->SetBranchAddress("reco4_x", &reco4_x , &b_reco4_x);
    tree->SetBranchAddress("res4_x", &res4_x , &b_res4_x);
    tree->SetBranchAddress("reco4_y", &reco4_y , &b_reco4_y);
    tree->SetBranchAddress("res4_y", &res4_y , &b_res4_y);
    tree->SetBranchAddress("z4", &z4 , &b_z4);
    tree->SetBranchAddress("dwc4_multiplicity", &averageHitMultiplicty4 , &b_averageHitMultiplicty4);


    tree->SetBranchAddress("XCET_021507_signal", &XCET_021507_signal, &b_XCET_021507_signal);
    tree->SetBranchAddress("XCET_021523_signal", &XCET_021523_signal, &b_XCET_021523_signal);
    tree->SetBranchAddress("N_scintillator_coincidence_timestamps", &scintillator_coincidences, &b_scintillator_coincidences);
    tree->SetBranchAddress("N_scintillator_veto_timestamps", &scintillator_vetos, &b_scintillator_vetos);
    tree->SetBranchAddress("valid_TS_MCP1", &valid_TS_MCP1, &b_valid_TS_MCP1);
    tree->SetBranchAddress("valid_TS_MCP2", &valid_TS_MCP2, &b_valid_TS_MCP2);
    tree->SetBranchAddress("TS_15PercentRise_MCP1", &TS_15PercentRise_MCP1, &b_TS_15PercentRise_MCP1);
    tree->SetBranchAddress("TS_15PercentRise_MCP2", &TS_15PercentRise_MCP2, &b_TS_15PercentRise_MCP2);
    tree->SetBranchAddress("TS_30PercentRise_MCP1", &TS_30PercentRise_MCP1, &b_TS_30PercentRise_MCP1);
    tree->SetBranchAddress("TS_30PercentRise_MCP2", &TS_30PercentRise_MCP2, &b_TS_30PercentRise_MCP2);
    tree->SetBranchAddress("TS_45PercentRise_MCP1", &TS_45PercentRise_MCP1, &b_TS_45PercentRise_MCP1);
    tree->SetBranchAddress("TS_45PercentRise_MCP2", &TS_45PercentRise_MCP2, &b_TS_45PercentRise_MCP2);
    tree->SetBranchAddress("TS_60PercentRise_MCP1", &TS_60PercentRise_MCP1, &b_TS_60PercentRise_MCP1);
    tree->SetBranchAddress("TS_60PercentRise_MCP2", &TS_60PercentRise_MCP2, &b_TS_60PercentRise_MCP2);
    tree->SetBranchAddress("amp_MCP1", &amp_MCP1, &b_amp_MCP1);
    tree->SetBranchAddress("amp_MCP2", &amp_MCP2, &b_amp_MCP2);

    tree->SetBranchAddress("TS_MCP1", &TS_MCP1, &b_TS_MCP1);
    tree->SetBranchAddress("TS_MCP2", &TS_MCP2, &b_TS_MCP2);
    tree->SetBranchAddress("TS_MCP1_to_last_falling_Edge", &TS_MCP1_to_last_falling_Edge, &b_TS_MCP1_to_last_falling_Edge);
    tree->SetBranchAddress("TS_MCP2_to_last_falling_Edge", &TS_MCP2_to_last_falling_Edge, &b_TS_MCP2_to_last_falling_Edge);

    loaded_run = -1;
}

void HGCalTBBeamWireChamberProducer::produce(edm::Event& event, const edm::EventSetup& iSetup) {

    //get the relevant event information
    edm::Handle<RunData> rd;
    event.getByToken(RunDataToken, rd);

    if (rd->run != loaded_run) {
        loadRun(rd->run);
    }

    std::unique_ptr<std::map<int, WireChamberData> > dwcs(new std::map<int, WireChamberData>);
    WireChamberData* dwc1 = new WireChamberData();
    WireChamberData* dwc2 = new WireChamberData();
    WireChamberData* dwc3 = new WireChamberData();
    WireChamberData* dwc4 = new WireChamberData();

    //set the RunData
    std::unique_ptr<RunData> rd_full(new RunData);

    rd_full->configuration = rd->configuration;
    rd_full->run = rd->run;
    rd_full->trigger = rd->trigger;
    rd_full->trigger_ts = rd->trigger_ts;
    rd_full->event = rd->event;
    rd_full->energy = rd->energy;
    rd_full->runType = rd->runType;
    rd_full->pdgID = rd->pdgID;

    if (rd->booleanUserRecords.has("hasDanger")) rd_full->booleanUserRecords.add("hasDanger", rd->booleanUserRecords.get("hasDanger"));

    if (rd->doubleUserRecords.has("trueEnergy")) rd_full->doubleUserRecords.add("trueEnergy", rd->doubleUserRecords.get("trueEnergy"));

    if (tree != NULL) {


        if (reco1_x_loaded.count(rd->event) == 0) {
            dwc1->goodMeasurement = dwc2->goodMeasurement = dwc3->goodMeasurement = dwc4->goodMeasurement = false;
            rd_full->booleanUserRecords.add("hasValidDWCMeasurement", false);
            rd_full->doubleUserRecords.add("triggerDeltaT_to_TDC", -999.);
        } else {
            //cross-check
            int event_nr = rd->event;

#ifdef DEBUG
            if (averageHitMultiplicty2_loaded[event_nr] > 100) {

                std::cout << "Run: " << rd->run << "   Event: " << event_nr << std::endl;
                std::cout << "GoodDWC: " << goodDWC_Measurement_loaded[event_nr] << "   triggerTimeDifference: " << triggerTimeDiff_loaded[event_nr] << std::endl;

                std::cout << "reco1_x: " << reco1_x_loaded[event_nr] << std::endl;
                std::cout << "reco1_y: " << reco1_y_loaded[event_nr] << std::endl;
                std::cout << "res1_x: " << res1_x_loaded[event_nr] << std::endl;
                std::cout << "res1_y: " << res1_y_loaded[event_nr] << std::endl;
                std::cout << "z1: " << z1_loaded[event_nr] << std::endl;
                std::cout << "multiplicity1: " << averageHitMultiplicty1_loaded[event_nr] << std::endl;

                std::cout << "reco2_x: " << reco2_x_loaded[event_nr] << std::endl;
                std::cout << "reco2_y: " << reco2_y_loaded[event_nr] << std::endl;
                std::cout << "res2_x: " << res2_x_loaded[event_nr] << std::endl;
                std::cout << "res2_y: " << res2_y_loaded[event_nr] << std::endl;
                std::cout << "z2: " << z2_loaded[event_nr] << std::endl;
                std::cout << "multiplicity2: " << averageHitMultiplicty2_loaded[event_nr] << std::endl;

                std::cout << "reco3_x: " << reco3_x_loaded[event_nr] << std::endl;
                std::cout << "reco3_y: " << reco3_y_loaded[event_nr] << std::endl;
                std::cout << "res3_x: " << res3_x_loaded[event_nr] << std::endl;
                std::cout << "res3_y: " << res3_y_loaded[event_nr] << std::endl;
                std::cout << "z3: " << z3_loaded[event_nr] << std::endl;
                std::cout << "multiplicity3: " << averageHitMultiplicty3_loaded[event_nr] << std::endl;

                std::cout << "reco4_x: " << reco4_x_loaded[event_nr] << std::endl;
                std::cout << "reco4_y: " << reco4_y_loaded[event_nr] << std::endl;
                std::cout << "res4_x: " << res4_x_loaded[event_nr] << std::endl;
                std::cout << "res4_y: " << res4_y_loaded[event_nr] << std::endl;
                std::cout << "z4: " << z4_loaded[event_nr] << std::endl << std::endl;
                std::cout << "multiplicity4: " << averageHitMultiplicty4_loaded[event_nr] << std::endl;
            }
#endif


            dwc1->ID = 1;
            dwc1->goodMeasurement_X = (reco1_x_loaded[event_nr] != -999);
            dwc1->goodMeasurement_Y = (reco1_y_loaded[event_nr] != -999);
            dwc1->goodMeasurement = (dwc1->goodMeasurement_X && dwc1->goodMeasurement_Y);
            dwc1->x = -reco1_x_loaded[event_nr];    //September 2017: Invert the sign in the x-coordinate
            dwc1->res_x = res1_x_loaded[event_nr];
            dwc1->y = reco1_y_loaded[event_nr];
            dwc1->res_y = res1_y_loaded[event_nr];
            dwc1->z = z1_loaded[event_nr];
            dwc1->averageHitMultiplicty = averageHitMultiplicty1_loaded[event_nr];


            dwc2->ID = 2;
            dwc2->goodMeasurement_X = (reco2_x_loaded[event_nr] != -999);
            dwc2->goodMeasurement_Y = (reco2_y_loaded[event_nr] != -999);
            dwc2->goodMeasurement = (dwc2->goodMeasurement_X && dwc2->goodMeasurement_Y);
            dwc2->x = -reco2_x_loaded[event_nr];    //September 2017: Invert the sign in the x-coordinate
            dwc2->res_x = res2_x_loaded[event_nr];
            dwc2->y = reco2_y_loaded[event_nr];
            dwc2->res_y = res2_y_loaded[event_nr];
            dwc2->z = z2_loaded[event_nr];
            dwc2->averageHitMultiplicty = averageHitMultiplicty2_loaded[event_nr];


            dwc3->ID = 1;
            dwc3->goodMeasurement_X = (reco3_x_loaded[event_nr] != -999);
            dwc3->goodMeasurement_Y = (reco3_y_loaded[event_nr] != -999);
            dwc3->goodMeasurement = (dwc3->goodMeasurement_X && dwc3->goodMeasurement_Y);
            dwc3->x = -reco3_x_loaded[event_nr];    //September 2017: Invert the sign in the x-coordinate
            dwc3->res_x = res3_x_loaded[event_nr];
            dwc3->y = reco3_y_loaded[event_nr];
            dwc3->res_y = res3_y_loaded[event_nr];
            dwc3->z = z3_loaded[event_nr];
            dwc3->averageHitMultiplicty = averageHitMultiplicty3_loaded[event_nr];


            dwc4->ID = 1;
            dwc4->goodMeasurement_X = (reco4_x_loaded[event_nr] != -999);
            dwc4->goodMeasurement_Y = (reco4_y_loaded[event_nr] != -999);
            dwc4->goodMeasurement = (dwc4->goodMeasurement_X && dwc4->goodMeasurement_Y);
            dwc4->x = -reco4_x_loaded[event_nr];    //September 2017: Invert the sign in the x-coordinate
            dwc4->res_x = res4_x_loaded[event_nr];
            dwc4->y = reco4_y_loaded[event_nr];
            dwc4->res_y = res4_y_loaded[event_nr];
            dwc4->z = z4_loaded[event_nr];
            dwc4->averageHitMultiplicty = averageHitMultiplicty4_loaded[event_nr];


            rd_full->booleanUserRecords.add("hasValidDWCMeasurement", (bool)goodDWC_Measurement_loaded[event_nr]);
            rd_full->doubleUserRecords.add("triggerDeltaT_to_TDC", triggerTimeDiff_loaded[event_nr]);


            //XCET
            rd_full->intUserRecords.add("XCET_021507_signal", XCET_021507_signal_loaded[event_nr]);
            rd_full->intUserRecords.add("XCET_021523_signal", XCET_021523_signal_loaded[event_nr]);
            //scintillators
            rd_full->intUserRecords.add("scintillator_coincidence_timestamps", scintillator_coincidences_loaded[event_nr]);
            rd_full->intUserRecords.add("scintillator_coincidence_timestamps", scintillator_vetos_loaded[event_nr]);
            //MCPs
            rd_full->intUserRecords.add("valid_TS_MCP1", valid_TS_MCP1_loaded[event_nr]);
            rd_full->intUserRecords.add("valid_TS_MCP2", valid_TS_MCP2_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_MCP1", TS_MCP1_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_MCP2", TS_MCP2_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_15PercentRise_MCP1", TS_15PercentRise_MCP1_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_15PercentRise_MCP2", TS_15PercentRise_MCP2_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_30PercentRise_MCP2", TS_30PercentRise_MCP2_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_30PercentRise_MCP1", TS_30PercentRise_MCP1_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_45PercentRise_MCP1", TS_45PercentRise_MCP1_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_45PercentRise_MCP2", TS_45PercentRise_MCP2_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_60PercentRise_MCP1", TS_60PercentRise_MCP1_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_60PercentRise_MCP2", TS_60PercentRise_MCP2_loaded[event_nr]);
            rd_full->doubleUserRecords.add("amp_MCP1", amp_MCP1_loaded[event_nr]);
            rd_full->doubleUserRecords.add("amp_MCP2", amp_MCP2_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_MCP1_to_last_falling_Edge", TS_MCP1_to_last_falling_Edge_loaded[event_nr]);
            rd_full->doubleUserRecords.add("TS_MCP2_to_last_falling_Edge", TS_MCP2_to_last_falling_Edge_loaded[event_nr]);
        }
    }
    (*dwcs)[0] = *dwc1;
    (*dwcs)[1] = *dwc2;
    (*dwcs)[2] = *dwc3;
    (*dwcs)[3] = *dwc4;

    event.put(std::move(dwcs), "DelayWireChambers");


    event.put(std::move(rd_full), "FullRunData");

}

void HGCalTBBeamWireChamberProducer::loadRun(int loading_run) {
#ifdef DEBUG
    std::cout << "Clearing run " << loaded_run << std::endl;
#endif
    reco1_x_loaded.clear(); reco1_y_loaded.clear(); z1_loaded.clear(); averageHitMultiplicty1_loaded.clear();
    reco2_x_loaded.clear(); reco2_y_loaded.clear(); z2_loaded.clear(); averageHitMultiplicty2_loaded.clear();
    reco3_x_loaded.clear(); reco3_y_loaded.clear(); z3_loaded.clear(); averageHitMultiplicty3_loaded.clear();
    reco4_x_loaded.clear(); reco4_y_loaded.clear(); z4_loaded.clear(); averageHitMultiplicty4_loaded.clear();
    XCET_021507_signal_loaded.clear(); XCET_021523_signal_loaded.clear(); scintillator_coincidences_loaded.clear(); scintillator_vetos_loaded.clear();
    valid_TS_MCP1_loaded.clear(); valid_TS_MCP2_loaded.clear();
    TS_15PercentRise_MCP1_loaded.clear(); TS_15PercentRise_MCP2_loaded.clear(); TS_30PercentRise_MCP2_loaded.clear(); TS_30PercentRise_MCP1_loaded.clear(); TS_45PercentRise_MCP1_loaded.clear(); TS_45PercentRise_MCP2_loaded.clear(); TS_60PercentRise_MCP1_loaded.clear(); TS_60PercentRise_MCP2_loaded.clear(); amp_MCP1_loaded.clear(); amp_MCP2_loaded.clear();
    TS_MCP1_loaded.clear(); TS_MCP2_loaded.clear(); TS_MCP1_to_last_falling_Edge_loaded.clear(); TS_MCP2_to_last_falling_Edge_loaded.clear();

    if (tree == NULL) return;

#ifdef DEBUG
    std::cout << "Loading run " << loading_run << std::endl;
#endif

    for (size_t i = 0; i < (size_t) tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (run != loading_run) continue;

#ifdef DEBUG
        std::cout << "Added event: " << eventId << std::endl;
#endif

        goodDWC_Measurement_loaded[eventId] = goodDWC_Measurement;
        triggerTimeDiff_loaded[eventId] = triggerTimeDiff;

        reco1_x_loaded[eventId] = reco1_x;
        res1_x_loaded[eventId] = res1_x;
        reco1_y_loaded[eventId] = reco1_y;
        res1_y_loaded[eventId] = res1_y;
        z1_loaded[eventId] = z1;
        averageHitMultiplicty1_loaded[eventId] = averageHitMultiplicty1;

        reco2_x_loaded[eventId] = reco2_x;
        res2_x_loaded[eventId] = res2_x;
        reco2_y_loaded[eventId] = reco2_y;
        res2_y_loaded[eventId] = res2_y;
        z2_loaded[eventId] = z2;
        averageHitMultiplicty2_loaded[eventId] = averageHitMultiplicty2;

        reco3_x_loaded[eventId] = reco3_x;
        res3_x_loaded[eventId] = res3_x;
        reco3_y_loaded[eventId] = reco3_y;
        res3_y_loaded[eventId] = res3_y;
        z3_loaded[eventId] = z3;
        averageHitMultiplicty3_loaded[eventId] = averageHitMultiplicty3;

        reco4_x_loaded[eventId] = reco4_x;
        res4_x_loaded[eventId] = res4_x;
        reco4_y_loaded[eventId] = reco4_y;
        res4_y_loaded[eventId] = res4_y;
        z4_loaded[eventId] = z4;
        averageHitMultiplicty4_loaded[eventId] = averageHitMultiplicty4;

        XCET_021507_signal_loaded[eventId] = XCET_021507_signal;
        XCET_021523_signal_loaded[eventId] = XCET_021523_signal;
        scintillator_coincidences_loaded[eventId] = scintillator_coincidences;
        scintillator_vetos_loaded[eventId] = scintillator_vetos;

        valid_TS_MCP1_loaded[eventId] = valid_TS_MCP1;
        valid_TS_MCP2_loaded[eventId] = valid_TS_MCP2;
        TS_MCP1_loaded[eventId] = TS_MCP1;
        TS_MCP2_loaded[eventId] = TS_MCP2;
        TS_15PercentRise_MCP1_loaded[eventId] = TS_15PercentRise_MCP1;
        TS_15PercentRise_MCP2_loaded[eventId] = TS_15PercentRise_MCP2;
        TS_30PercentRise_MCP2_loaded[eventId] = TS_30PercentRise_MCP2;
        TS_30PercentRise_MCP1_loaded[eventId] = TS_30PercentRise_MCP1;
        TS_45PercentRise_MCP1_loaded[eventId] = TS_45PercentRise_MCP1;
        TS_45PercentRise_MCP2_loaded[eventId] = TS_45PercentRise_MCP2;
        TS_60PercentRise_MCP1_loaded[eventId] = TS_60PercentRise_MCP1;
        TS_60PercentRise_MCP2_loaded[eventId] = TS_60PercentRise_MCP2;
        amp_MCP1_loaded[eventId] = amp_MCP1;
        amp_MCP2_loaded[eventId] = amp_MCP2;
        TS_MCP1_to_last_falling_Edge_loaded[eventId] = TS_MCP1_to_last_falling_Edge;
        TS_MCP2_to_last_falling_Edge_loaded[eventId] = TS_MCP2_to_last_falling_Edge;
    }
#ifdef DEBUG
    std::cout << "Loaded run " << loading_run << std::endl;
#endif

    loaded_run = loading_run;
}

void HGCalTBBeamWireChamberProducer::endJob() {
    if (tree != NULL) delete tree;
    if (rootFile != NULL) delete rootFile;

    reco1_x_loaded.clear(); reco1_y_loaded.clear(); z1_loaded.clear(); averageHitMultiplicty1_loaded.clear();
    reco2_x_loaded.clear(); reco2_y_loaded.clear(); z2_loaded.clear(); averageHitMultiplicty1_loaded.clear();
    reco3_x_loaded.clear(); reco3_y_loaded.clear(); z3_loaded.clear(); averageHitMultiplicty1_loaded.clear();
    reco4_x_loaded.clear(); reco4_y_loaded.clear(); z4_loaded.clear(); averageHitMultiplicty1_loaded.clear();

    res1_x_loaded.clear(); res1_y_loaded.clear();
    res2_x_loaded.clear(); res2_y_loaded.clear();
    res3_x_loaded.clear(); res3_y_loaded.clear();
    res4_x_loaded.clear(); res4_y_loaded.clear();
}

DEFINE_FWK_MODULE(HGCalTBBeamWireChamberProducer);
