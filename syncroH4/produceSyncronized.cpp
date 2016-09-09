//g++ -Wall -o produceSyncronized `root-config --cflags --glibs` H4treeWC.cc H4treeRecoWC.cc produceSyncronized.cpp  
// ./produceSyncronized   H4run pathToH4_run
//H4run = folder containing the nSpill.root

#include "TChain.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"

#include "H4treeWC.h"
#include "H4treeRecoWC.h"



int main(int argc, char** argv){

  std::string H4_run = argv[1];
  std::string pathToH4_run = argv[2];

  std::cout << " reading H4DAQ run = " << H4_run << std::endl;

  TChain* tree = new TChain("H4tree", ("H4treeRecoWCOut_"+H4_run+".root").c_str());
  if(argc < 3) tree->Add((H4_run+"/*.root").c_str());
  else tree->Add((pathToH4_run+"/"+H4_run+"/*.root").c_str());


  H4treeRecoWC* H4Recotree = new H4treeRecoWC(tree);

  int nEntries_beam = H4Recotree->TotEntries();
  std::cout << " >> entries = " << nEntries_beam << std::endl;

  
  ULong64_t evtTime_prev = 0;
  ULong64_t evtTime_deltaPrev = 0;

  int countInterSpill = 0;
  int totInterSpill = 0;

  std::vector<double> evtPerSpill;
  std::vector<double> dT;
  std::vector<double> WC_x0;
  std::vector<double> WC_x1;
  std::vector<double> WC_y0;
  std::vector<double> WC_y1;

  evtPerSpill.clear();
  
  for(int iBeam=0; iBeam<nEntries_beam; ++iBeam){

    if(iBeam == 0) evtTime_deltaPrev = 0;
    else  evtTime_prev = H4Recotree->evtTime[0];
    
    H4Recotree->GetEntry(iBeam);

    if(iBeam != 0) evtTime_deltaPrev = H4Recotree->evtTime[0] - evtTime_prev;
    H4Recotree->FillTDC();

    // std::cout << H4Recotree->evtNumber << " \t " << evtTime_deltaPrev << " \t " 
    // 	      << H4Recotree->wc_recox_[0] << " \t " << H4Recotree->wc_recoy_[0] << " \t " 
    // 	      << H4Recotree->wc_recox_[1] << " \t " << H4Recotree->wc_recoy_[1] << std::endl;


    //    if(iBeam > 0)
    //    std::cout << iBeam << " \t " << evtTime_deltaPrev << std::endl;

    H4Recotree->FillEvent();
    dT.push_back(evtTime_deltaPrev);
    WC_x0.push_back(H4Recotree->wc_recox_[0]);
    WC_x1.push_back(H4Recotree->wc_recox_[1]);
    WC_y0.push_back(H4Recotree->wc_recoy_[0]);
    WC_y1.push_back(H4Recotree->wc_recoy_[1]);

    if(evtTime_deltaPrev > 1.e+06){
      evtPerSpill.push_back(countInterSpill);
      countInterSpill = 0;
    }
    if(iBeam == nEntries_beam - 1){
      evtPerSpill.push_back(countInterSpill);
    }

    ++countInterSpill;
    ++totInterSpill;
    
  } //raw ntuple with WC coded info


  //  return 100 ;

  delete H4Recotree;

  unsigned int nSpill = 0;
  unsigned int evtSpill = evtPerSpill.at(nSpill);

  // write out txt
  std::ofstream outFile(("WC_H4Run"+H4_run+".txt").c_str(), std::ios::out);
  outFile << "H4run " << atoi(H4_run.c_str()) << " evtSpill " << evtSpill << std::endl;

  
  for(unsigned int iEv = 0; iEv < dT.size(); ++iEv){

    if(iEv  == evtSpill && iEv+1 != dT.size()){
      outFile << "H4run " << atoi(H4_run.c_str()) << " evtSpill " << evtPerSpill.at(nSpill) << std::endl;
      //std::cout << " fine spill " << nSpill << " nEvt = " << evtSpill << " dEvt = " << evtPerSpill.at(nSpill) << " iEv = " << iEv << " dT = " << dT.at(iEv) << std::endl;
      if(nSpill < evtPerSpill.size()-1){
	++nSpill;
	evtSpill += evtPerSpill.at(nSpill);
      }
    }  
    outFile << WC_x0.at(iEv) << " \t " << WC_x1.at(iEv) << " \t " << WC_y0.at(iEv) << " \t " << WC_y1.at(iEv) << std::endl;
  }
  
  outFile.close();


  std::cout << " END " << std::endl;
}
