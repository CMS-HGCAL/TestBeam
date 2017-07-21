#ifndef HGCALTBRECHITPRODUCER_H
#define HGCALTBRECHITPRODUCER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include <iostream>
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBWireChamberData.h"

#include "TFile.h"
#include "TTree.h"

class HGCalTBBeamWireChamberProducer : public edm::EDProducer{

  public:
    HGCalTBBeamWireChamberProducer(const edm::ParameterSet&);
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
 
  private:
    edm::EDGetTokenT<RunData> RunDataToken; 
    edm::EDGetTokenT<WireChambers> MWCToken;

    std::string inputFile;
    TFile* rootFile;
    TTree* tree;


    int run, eventId, goodDWC_Measurement;
    double reco1_x, reco1_y, z1, reco2_x, reco2_y, z2, reco3_x, reco3_y, z3, reco4_x, reco4_y, z4;    
    TBranch                   *b_run;   
    TBranch                   *b_event;   
    TBranch                   *b_goodDWC_Measurement;   
    TBranch                   *b_reco1_x;   
    TBranch                   *b_reco1_y;   
    TBranch                   *b_z1;   
    TBranch                   *b_reco2_x;   
    TBranch                   *b_reco2_y;   
    TBranch                   *b_z2;   
    TBranch                   *b_reco3_x;   
    TBranch                   *b_reco3_y;   
    TBranch                   *b_z3;   
    TBranch                   *b_reco4_x;   
    TBranch                   *b_reco4_y;   
    TBranch                   *b_z4;       


    int loaded_run;
    std::map<int, double> reco1_x_loaded;
    std::map<int, double> reco1_y_loaded;
    std::map<int, double> z1_loaded;

    std::map<int, double> reco2_x_loaded;
    std::map<int, double> reco2_y_loaded;
    std::map<int, double> z2_loaded;
    
    std::map<int, double> reco3_x_loaded;
    std::map<int, double> reco3_y_loaded;
    std::map<int, double> z3_loaded;

    std::map<int, double> reco4_x_loaded;
    std::map<int, double> reco4_y_loaded;
    std::map<int, double> z4_loaded;

    void loadRun(int loading_run);
};
#endif
