//author: Thorben Quast
//date: 18th May 2018
//producer reading reconstructed DATURA telescope data and additing to the root file
//v1 (18th May 2018) also implements a preliminary straight line tracking to compute reference points on the sensors

#ifndef HGCALTBDATURATELESCOPEPRODUCER_H
#define HGCALTBDATURATELESCOPEPRODUCER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include <iostream>
#include "HGCal/DataFormats/interface/HGCalTBRunData.h"
#include "HGCal/DataFormats/interface/HGCalTBDATURATelescopeData.h"

//for sttraight line tracking
#include "HGCal/Reco/interface/PositionResolutionHelpers.h"

#include "TFile.h"
#include "TTree.h"

class HGCalTBDATURATelescopeProducer : public edm::EDProducer{

  public:
    HGCalTBDATURATelescopeProducer(const edm::ParameterSet&);
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
 
  private:
    edm::EDGetTokenT<RunData> RunDataToken; 
    
    std::string outputCollectionName;

    std::string inputFile;
    TFile* rootFile;
    TTree* tree;
    int SkipFirstNEvents;

     
    TBranch                   *b_event;   
    TBranch                   *b_Ntracks;  
    TBranch                   *b_chi2;  
    std::map<int, TBranch*> b_clusterX;
    std::map<int, TBranch*> b_clusterY;
    std::map<int, TBranch*> b_clusterZ;
    std::map<int, TBranch*> b_absorber;


    int tree_eventID;
    int tree_Ntracks;
    std::vector<double> *tree_track_chi2;
    std::map<int, std::vector<double>* > tree_clusterX;       //in mm
    std::map<int, std::vector<double>* > tree_clusterY;       //in mm
    std::map<int, std::vector<double>* > tree_clusterZ;       //in mm
    std::map<int, std::vector<double>* > tree_absorber;       //in X0

    //for computation of reference positions at each layer
    std::string m_layerPositionFile;
    std::map<int, double> layerPositions;
    
};
#endif
