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
#include "../MCanalysis/HGCalTBPlots.C"

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
  //  float MIP2GeV_sim = EtoMip * 63.28e-06 / 55.16e-06;  // recalibrate to mean 500MeV muon response
  float MIP2GeV_sim = EtoMip * 51.91e-06 / 55.16e-06;  // recalibrate to mean 500MeV muon response

  float G4Escale = 1.; // 0.83;   // >>> fix scale in data mc for Si


  TChain* t = new TChain("HGCalTBAnalyzer/HGCTB");
  
  if(nLayer == 8 && config == 1){
    std::string nameFolder = "";
    if(energyEle == 20) nameFolder = "161116_234146";
    if(energyEle == 32) nameFolder = "161117_101400";
    if(energyEle == 70) nameFolder = "161117_082845";
    if(energyEle == 100) nameFolder = "161116_232458";
    #if(energyEle == 150) nameFolder = "161116_232510";
    if(energyEle == 200) nameFolder = "161116_234226";
    if(energyEle == 250) nameFolder = "161116_234238";
    
 for(int iFo=0; iFo<2; ++iFo){
 for(int iF=0; iF<1000; ++iF){
 t->Add(Form(("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/8moduleIv4_PhysicsList_FTFP_BERT_EMM_beamposition_spread_Momspread_realistic_nosecKapton/mc/CRAB_PrivateMC/crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
 }
 }
  }
  else if(nLayer == 8 && config == 2){
    std::string nameFolder = "";
    if(energyEle == 20) nameFolder = "161119_013149";
    if(energyEle == 32) nameFolder = "161119_013155";
    if(energyEle == 70) nameFolder = "161119_013201";
    if(energyEle == 100) nameFolder = "161119_013207";
    //if(energyEle == 150) nameFolder = "";
    if(energyEle == 200) nameFolder = "161119_013217";
    if(energyEle == 250) nameFolder = "161119_013223";

 for(int iFo=0; iFo<2; ++iFo){
 for(int iF=0; iF<1000; ++iF){
 t->Add(Form(("root://eoscms.cern.ch//store/group/upgrade/HGCAL/simulation/8moduleIIv5_PhysicsList_FTFP_BERT_EMM_beamposition_spread_Momspread_realistic_nosecKapton/mc/CRAB_PrivateMC//crab_Ele%dGeV/"+nameFolder+"/000%d/TBGenSim_%d.root").c_str(), energyEle, iFo, iF));
 }
 }    
  }
  else{
    std::string nameFolder = "";
    if(energyEle == 20) nameFolder = "161031_180325";
    //if(energyEle == 35) nameFolder = "161031_180332";
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

  TH1F* enLayer_e19[28];
  TH1F* enLayerCorr_e19[28];
  TH1F* enLayerMip_e19[28];
  TH1F* enLayer_eAll[28];
  TH1F* enLayerCorr_eAll[28];
  TH1F* enLayerMip_eAll[28];
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer_e19[i] = new TH1F(Form("enLayer%d_e19", i+1), "", 5000, 0., 0.5);
    enLayerCorr_e19[i] = new TH1F(Form("enLayerCorr%d_e19", i+1), "", 5000, 0., 10.);
    enLayerMip_e19[i] = new TH1F(Form("enLayerMip%d_e19", i+1), "", 10000, 0., 10000.);
    enLayer_eAll[i] = new TH1F(Form("enLayer%d_eAll", i+1), "", 5000, 0., 0.5);
    enLayerCorr_eAll[i] = new TH1F(Form("enLayerCorr%d_eAll", i+1), "", 5000, 0., 10.);
    enLayerMip_eAll[i] = new TH1F(Form("enLayerMip%d_eAll", i+1), "", 10000, 0., 10000.);
  }

  TH1F* recoEnergy_eAll = new TH1F("recoEnergy_eAll", "", 5000, 0., 500.);
  TH1F* recoEnergyRel_eAll = new TH1F("recoEnergyRel_eAll", "", 1000, 0., 2.);
  TGraphErrors* enLayer_MIP_eAll = new TGraphErrors();

  TH1F* recoEnergy_e19 = new TH1F("recoEnergy_e19", "", 5000, 0., 500.);
  TH1F* recoEnergyRel_e19 = new TH1F("recoEnergyRel_e19", "", 1000, 0., 2.);
  TGraphErrors* enLayer_MIP_e19 = new TGraphErrors();

  double xBeam, yBeam;
  std::vector<float>* simHitLayEn2E = 0;
  std::vector<float>* simHitCellEnE = 0;
  std::vector<unsigned int>* simHitCellIdE = 0;
  t->SetBranchAddress("simHitLayEn2E", &simHitLayEn2E);
  t->SetBranchAddress("simHitCellEnE", &simHitCellEnE);
  t->SetBranchAddress("simHitCellIdE", &simHitCellIdE);
  t->SetBranchAddress("xBeam", &xBeam);
  t->SetBranchAddress("yBeam", &yBeam);

  HexTopology ht1 = false;

  //looping over entries
  for(int iE=0; iE<totEntries; ++iE){

    t->GetEntry(iE);
    //    std::cout << " simHitLayEn2E->size() = " << simHitLayEn2E->size() << std::endl;
    if(abs(xBeam) >= 2 || abs(yBeam) >= 2) continue;
    float totEnergy_e19 = 0.;
    float totEnergy_eAll = 0.;
    /*
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
    */

    //////// check single cells
    std::vector<unsigned int> CellId;
    std::vector<float> CellE;
    for(int ii=0; ii<simHitCellIdE->size(); ++ii){
      CellId.push_back(simHitCellIdE->at(ii));
      CellE.push_back(simHitCellEnE->at(ii));
    }

    for(int iS=0; iS<simHitLayEn2E->size(); ++iS){   
      unsigned int locMaxId = ht1.localMax( (CellId), (CellE), iS+1);
      double clusterE19 = ht1.cluster((CellId), (CellE), locMaxId, 2, MIP2GeV_sim, 2);
      float localE_e19 = clusterE19 / G4Escale * ( 1./MIP2GeV_sim * weights2GeV * weights2MIP* dEdX_weights[iS] + 1.);                                                     
      
      double allcell = ht1.cluster( (CellId), (CellE), locMaxId, 7, MIP2GeV_sim, 2);
      float localE_eAll = allcell / G4Escale * ( 1./MIP2GeV_sim * weights2GeV * weights2MIP* dEdX_weights[iS] + 1.);                                                     
      
      enLayerCorr_e19[iS]->Fill(localE_e19);                                                                                                                                             
      enLayerMip_e19[iS]->Fill(clusterE19 / G4Escale * 1./MIP2GeV_sim);                                                                                                    
      totEnergy_e19 += localE_e19; 

      enLayerCorr_eAll[iS]->Fill(localE_eAll);                                                                                                                                             
      enLayerMip_eAll[iS]->Fill(allcell / G4Escale * 1./MIP2GeV_sim);                                                                                                    
      totEnergy_eAll += localE_eAll; 
    }
    recoEnergy_e19->Fill(totEnergy_e19);                   
    recoEnergyRel_e19->Fill(totEnergy_e19/(1.*energyEle)); 
    recoEnergy_eAll->Fill(totEnergy_eAll);                   
    recoEnergyRel_eAll->Fill(totEnergy_eAll/(1.*energyEle)); 
  }


  float X0value = 0;
  for(int i=0; i<nLayer; ++i){
    X0value += X0val[i];
    enLayer_MIP_eAll->SetPoint(i+1, X0value, enLayerMip_eAll[i]->GetMean());
    enLayer_MIP_eAll->SetPointError(i+1, 0, enLayerMip_eAll[i]->GetRMS()/enLayerMip_eAll[i]->GetEntries());
    enLayer_MIP_e19->SetPoint(i+1, X0value, enLayerMip_e19[i]->GetMean());
    enLayer_MIP_e19->SetPointError(i+1, 0, enLayerMip_e19[i]->GetRMS()/enLayerMip_e19[i]->GetEntries());
  }

  TFile* outF;
  if(nLayer == 28) outF = new TFile(Form("energyLayer_ROOT/outF_28layers_Ele%d.root", energyEle), "recreate");
  else if(config == 1) outF = new TFile(Form("energyLayer_ROOT/outF_Ele%d.root", energyEle), "recreate");
  else if(config == 2) outF = new TFile(Form("energyLayer_ROOT/outF_endShower_Ele%d.root", energyEle), "recreate");
  outF->cd();
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer_eAll[i]->Write();
    enLayerCorr_eAll[i]->Write();
    enLayerMip_eAll[i]->Write();

    enLayer_e19[i]->Write();
    enLayerCorr_e19[i]->Write();
    enLayerMip_e19[i]->Write();
  }
  enLayer_MIP_eAll->Write("enLayer_MIP_eAll");
  recoEnergy_eAll->Write();
  recoEnergyRel_eAll->Write();

  enLayer_MIP_e19->Write("enLayer_MIP_e19");
  recoEnergy_e19->Write();
  recoEnergyRel_e19->Write();
  outF->Close();
}
