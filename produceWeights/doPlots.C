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

void doPlots(int nLayer, int energyEle, int config){
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

  // pion MPV = 55.16;
  // muon MPV = 51.91;
  // muon Mean = 63.28;

  float EtoMip = 55.16e-06;  // 125pion MPV response

  float weights2GeV = 1.e-03;
  float weights2MIP = 1.;         // weights with MIP from PDG
  float MIP2GeV_sim = EtoMip * 63.28e-06 / 55.16e-06;  // recalibrate to mean 500MeV muon response

  float G4Escale = 1.; // 0.83;   // >>> fix scale in data mc for Si


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
  else{
    std::string nameFolder = "";
    if(energyEle == 20) nameFolder = "161031_180325";
    if(energyEle == 35) nameFolder = "161031_180332";
    if(energyEle == 70) nameFolder = "161031_180339";
    if(energyEle == 100) nameFolder = "161031_180346";
    if(energyEle == 150) nameFolder = "161031_180353";
    if(energyEle == 200) nameFolder = "161031_180401";
    if(energyEle == 250) nameFolder = "161101_095433";


for(int iFo=0; iFo<2; ++iFo){
for(int iF=0; iF<1000; ++iF){
t->Add(Form(("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/28moduleIv3_PhysicsList_FTFP_BERT_EMM/mc/CRAB_PrivateMC/crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
 }
 }
  }

  double totEntries = t->GetEntries();
  std::cout << " tree read => entries = " << totEntries << std::endl;

  TH1F* enLayer[28];
  TH1F* enLayerCorr[28];
  TH1F* enLayerMip[28];
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer[i] = new TH1F(Form("enLayer%d", i+1), "", 5000, 0., 0.5);
    enLayerCorr[i] = new TH1F(Form("enLayerCorr%d", i+1), "", 5000, 0., 10.);
    enLayerMip[i] = new TH1F(Form("enLayerMip%d", i+1), "", 10000, 0., 10000.);
  }

  TH1F* recoEnergy = new TH1F("recoEnergy", "", 5000, 0., 500.);
  TH1F* recoEnergyRel = new TH1F("recoEnergyRel", "", 1000, 0., 2.);
  TGraphErrors* enLayer_MIP = new TGraphErrors();

  double xBeam, yBeam;
  std::vector<float>* simHitLayEn2E = 0;
  t->SetBranchAddress("simHitLayEn2E", &simHitLayEn2E);
  t->SetBranchAddress("xBeam", &xBeam);
  t->SetBranchAddress("yBeam", &yBeam);


  //looping over entries
  for(int iE=0; iE<totEntries; ++iE){

    t->GetEntry(iE);
    //    std::cout << " simHitLayEn2E->size() = " << simHitLayEn2E->size() << std::endl;
    if(abs(xBeam) >= 2 || abs(yBeam) >= 2) continue;
    float totEnergy = 0.;
    for(int iS=0; iS<simHitLayEn2E->size(); ++iS){
      enLayer[iS]->Fill(simHitLayEn2E->at(iS));
      
      //      std::cout << " E = " simHitLayEn2E->at(iS) << << std::endl;
      float localE = simHitLayEn2E->at(iS) / G4Escale * ( 1./MIP2GeV_sim * weights2GeV * weights2MIP* dEdX_weights[iS] + 1.);

      enLayerCorr[iS]->Fill(localE);
      enLayerMip[iS]->Fill(simHitLayEn2E->at(iS) / G4Escale * 1./MIP2GeV_sim);
      totEnergy += localE;
    }
    recoEnergy->Fill(totEnergy);
    recoEnergyRel->Fill(totEnergy/(1.*energyEle));
  }


  float X0value = 0;
  for(int i=0; i<nLayer; ++i){
    X0value += X0val[i];
    enLayer_MIP->SetPoint(i+1, X0value, enLayerMip[i]->GetMean());
    enLayer_MIP->SetPointError(i+1, 0, enLayerMip[i]->GetRMS()/enLayerMip[i]->GetEntries());
  }

  TFile* outF;
  if(nLayer == 28) outF = new TFile(Form("energyLayer_ROOT/outF_28layers_Ele%d.root", energyEle), "recreate");
  else if(config == 1) outF = new TFile(Form("energyLayer_ROOT/outF_Ele%d.root", energyEle), "recreate");
  else if(config == 2) outF = new TFile(Form("energyLayer_ROOT/outF_endShower_Ele%d.root", energyEle), "recreate");
  outF->cd();
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer[i]->Write();
    enLayerCorr[i]->Write();
    enLayerMip[i]->Write();
  }
  enLayer_MIP->Write("enLayer_MIP");
  recoEnergy->Write();
  recoEnergyRel->Write();
  outF->Close();


}
