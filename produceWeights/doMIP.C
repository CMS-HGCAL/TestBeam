#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"

#include "/afs/cern.ch/work/a/amartell/ES_P5/CMSSW_7_4_6_patch6/src/MyAnalyzer/ESAnalyzer/macroAnalysis/langausfit.h"

void doMIP(int nLayer,  int config, int energyEle=0){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptTitle(0);

  float dEdX_weights[28];
  float X0val[28];

  if(nLayer == 28){
    float dEdX_weights_28[28] = {7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 7.203, 8.087, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 9.27, 10.817, 12.789, 12.789, 12.789, 12.789, 12.789, 12.789, 12.789, 8.148};
    float X0val_28[28] = {0.574, 0.699, 0.574, 0.699, 0.574, 0.699, 0.574, 0.699, 0.574, 0.699, 0.802, 0.977, 0.802, 0.977, 0.802, 0.977, 0.802, 0.977, 0.802, 0.977, 1.202, 1.439, 1.202, 1.439, 1.202, 1.439, 1.202, 1.439};

    for(int j=0; j<28; ++j){
      dEdX_weights[j] = dEdX_weights_28[j];
      X0val[j] = X0val_28[j];
    }
  }
  else{
    if(config == 1){
      //    float dEdX_weights_8[28] = {33.151, 13.261, 14.247, 9.885, 9.844, 9.844, 16.416, 28.296, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      float dEdX_weights_8[28] = {33.074, 13.184, 14.17, 9.788, 9.766, 9.766, 16.339, 14.129, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      float X0val_8[28] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      
      for(int j=0; j<28; ++j){
	dEdX_weights[j] = dEdX_weights_8[j];
	X0val[j] = X0val_8[j];
      }
    }
    else if(config == 2){
     float dEdX_weights_8[28] = {35.866, 30.864, 28.803, 23.095, 20.657, 19.804, 36.322, 27.451, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     float X0val_8[28] = {5.048, 3.412, 3.412, 2.866, 2.512, 1.625, 2.368, 6.021, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      
      for(int j=0; j<28; ++j){
	dEdX_weights[j] = dEdX_weights_8[j];
	X0val[j] = X0val_8[j];
      }
    }
  }


  TChain* t = new TChain("HGCalTBAnalyzer/HGCTB");
  
  if(nLayer == 8 && config == 1){
    std::string nameFolder = "";
    if(energyEle == 20) nameFolder = "161101_094655";
    if(energyEle == 35) nameFolder = "161031_180636";
    if(energyEle == 70) nameFolder = "161031_180644";
    if(energyEle == 100) nameFolder = "161101_094708";
    if(energyEle == 200) nameFolder = "161031_180705";
    if(energyEle == 250) nameFolder = "161031_180712";
    
 for(int iFo=0; iFo<2; ++iFo){
 for(int iF=0; iF<1000; ++iF){
 t->Add(Form(("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/8moduleIv3_PhysicsList_FTFP_BERT_EMM/mc/CRAB_PrivateMC/crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
 }
 }
  }
  else if(nLayer == 8 && config == 2){
    std::string nameFolder = "";
    if(energyEle == 20) nameFolder = "161108_231750";
    if(energyEle == 32) nameFolder = "161108_231756";
    if(energyEle == 70) nameFolder = "161109_074724";
    if(energyEle == 100) nameFolder = "161108_231811";
    if(energyEle == 200) nameFolder = "161109_130453";
    if(energyEle == 250) nameFolder = "161108_231833";

 for(int iFo=0; iFo<2; ++iFo){
 for(int iF=0; iF<1000; ++iF){
 t->Add(Form(("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/8moduleIIv4_PhysicsList_FTFP_BERT_EMM/mc/CRAB_PrivateMC/crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
 }
 }    
  }
  else if(nLayer == 28 && config == 1){
    for(int iFo=0; iFo<2; ++iFo){
      for(int iF=0; iF<1000; ++iF){
	//	t->Add(Form("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/28moduleIv3_PhysicsList_FTFP_BERT_EMM/mc/CRAB_PrivateMC/crab_Muon500MeV/161101_141550/000%d/TBGenSim_%d.root", iFo, iF));
t->Add(Form("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/28moduleIv3_PhysicsLst_QGSP_FTFP_BERT/mc/CRAB_PrivateMC/crab_Muon500MeV/161101_141550/000%d/TBGenSim_%d.root", iFo, iF));
      }
    }
  }
  else if(nLayer == 28 && config == 2){
    for(int iFo=0; iFo<2; ++iFo){
      for(int iF=0; iF<1000; ++iF){
	t->Add(Form("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/28moduleIv3_PhysicsList_FTFP_BERT_EMM/mc/CRAB_PrivateMC/crab_Pi125GeV/161101_095441/000%d/TBGenSim_%d.root", iFo, iF));
      }
    }
  }

  double totEntries = t->GetEntries();
  std::cout << " tree read => entries = " << totEntries << std::endl;
  
  TH1F* enLayer[28];
  
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer[i] = new TH1F(Form("enLayer%d", i+1), "", 5000, 0., 150.e-06);

  }


  TGraphErrors* enLayer_Mean = new TGraphErrors();
  TGraphErrors* enLayer_MPV = new TGraphErrors();

  double xBeam, yBeam;
  std::vector<float>* simHitLayEn2E = 0;
  t->SetBranchAddress("simHitLayEn2E", &simHitLayEn2E);
  t->SetBranchAddress("xBeam", &xBeam);
  t->SetBranchAddress("yBeam", &yBeam);


  //looping over entries
  for(int iE=0; iE<totEntries; ++iE){

    t->GetEntry(iE);
    nLayer = simHitLayEn2E->size();
    //    std::cout << " simHitLayEn2E->size() = " << simHitLayEn2E->size() << std::endl;
    if(abs(xBeam) >= 2 || abs(yBeam) >= 2) continue;
    float totEnergy = 0.;
    for(int iS=0; iS<simHitLayEn2E->size(); ++iS){
      if(simHitLayEn2E->at(iS) > 0)      enLayer[iS]->Fill(simHitLayEn2E->at(iS));
    }
  }

  TF1* hL = new TF1("hL", "landau", 0, 100);
  //  float X0value = 0;
  for(int i=0; i<nLayer; ++i){
    enLayer[i]->Fit("hL");

    //    X0value += X0val[i];
    enLayer_MPV->SetPoint(i+1, i+1, hL->GetParameter(1));
    enLayer_MPV->SetPointError(i+1, 0, hL->GetParError(1));

    enLayer_Mean->SetPoint(i+1, i+1, enLayer[i]->GetMean());
    enLayer_Mean->SetPointError(i+1, 0, enLayer[i]->GetRMS()/enLayer[i]->GetEntries());
  }

  TFile* outF;
  if(config == 1) outF = new TFile("energyLayer_ROOT/outF_MIPmuon_28layers.root", "recreate");
  if(config == 2) outF = new TFile("energyLayer_ROOT/outF_MIPpion_28layers.root", "recreate");
  /*
  else if(nLayer == 8 && config == 1) outF = new TFile(Form("energyLayer_ROOT/outF_Ele%d.root", energyEle), "recreate");
  else if(nLayer == 8 && config == 2) outF = new TFile(Form("energyLayer_ROOT/outF_endShower_Ele%d.root", energyEle), "recreate");
  */
  outF->cd();
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer[i]->Write();
  }
  enLayer_Mean->Write("enLayer_Mean");
  enLayer_MPV->Write("enLayer_MPV");
  outF->Close();


}
