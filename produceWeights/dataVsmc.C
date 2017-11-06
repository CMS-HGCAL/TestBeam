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



void dataVsmc(){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

               
  int nLayers = 8;
  int config = 2;
  int iColors[7] = {kRed, kBlue, kGreen+2, kRed, kBlue, kRed, kBlue};
  //  int energyP[3] = {20, 70, 100, 200, 250};
  
  TFile* inF[4];
  inF[0] = TFile::Open("SIM_resolution_layers8_config1.root");
  inF[1] = TFile::Open("SIM_resolution_layers8_config2.root");
  inF[2] = TFile::Open("DATA_resolution_layers8_config1.root");
  inF[3] = TFile::Open("DATA_resolution_layers8_config2.root");


  TGraphErrors* tgCfg1[7];
  TGraphErrors* tgCfg2[7];

  tgCfg1[0] = (TGraphErrors*)inF[0]->Get("resolution_GeV");
  tgCfg1[1] = (TGraphErrors*)inF[2]->Get("resolution_GeV");
  tgCfg1[2] = (TGraphErrors*)inF[2]->Get("resolution_Mip");
  tgCfg1[3] = (TGraphErrors*)inF[0]->Get("linerity_GeV");
  tgCfg1[4] = (TGraphErrors*)inF[2]->Get("linearity_GeV");
  tgCfg1[5] = (TGraphErrors*)inF[0]->Get("mean_GeV");
  tgCfg1[6] = (TGraphErrors*)inF[2]->Get("mean_GeV");

  tgCfg2[0] = (TGraphErrors*)inF[1]->Get("resolution_GeV");
  tgCfg2[1] = (TGraphErrors*)inF[3]->Get("resolution_GeV");
  tgCfg2[2] = (TGraphErrors*)inF[3]->Get("resolution_Mip");
  tgCfg2[3] = (TGraphErrors*)inF[1]->Get("linerity_GeV");
  tgCfg2[4] = (TGraphErrors*)inF[3]->Get("linearity_GeV");
  tgCfg2[5] = (TGraphErrors*)inF[1]->Get("mean_GeV");
  tgCfg2[6] = (TGraphErrors*)inF[3]->Get("mean_GeV");


  for(int iT=0; iT<7; ++iT){
    tgCfg1[iT]->SetMarkerStyle(20);
    tgCfg2[iT]->SetMarkerStyle(20);

    tgCfg1[iT]->SetMarkerColor(iColors[iT]);
    tgCfg2[iT]->SetMarkerColor(iColors[iT]);
  }

  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.45,0.35,0.65,0.55,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->SetHeader("6 X_{0} - 15 X_{0}");
  legTGM->AddEntry(tgCfg1[0], "sim GeV", "pl");
  legTGM->AddEntry(tgCfg1[1], "data GeV", "pl");
  //  legTGM->AddEntry(tgCfg1[2], "data Mip", "l");

  TLegend *legTGM2 = new TLegend(0.45,0.35,0.65,0.55,NULL,"brNDC");
  legTGM2->SetTextFont(42);
  legTGM2->SetFillColor(kWhite);
  legTGM2->SetLineColor(kWhite);
  legTGM2->SetShadowColor(kWhite);
  legTGM2->SetHeader("5 X_{0} - 27 X_{0}");
  legTGM2->AddEntry(tgCfg1[0], "sim GeV", "pl");
  legTGM2->AddEntry(tgCfg1[1], "data GeV", "pl");
  //  legTGM->AddEntry(tgCfg1[2], "data Mip", "l");

  

  TLegend *leg = new TLegend(0.45,0.65,0.85,0.85,NULL,"brNDC");
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->SetHeader("6 X_{0} - 15 X_{0}");
  leg->AddEntry(tgCfg1[1], "data GeV", "pl");
  leg->AddEntry(tgCfg1[2], "data Mip (no weights)", "pl");   


  TLegend *leg2 = new TLegend(0.45,0.65,0.85,0.85,NULL,"brNDC");
  leg2->SetTextFont(42);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->SetShadowColor(kWhite);
  leg2->SetHeader("5 X_{0} - 27 X_{0}");
  leg2->AddEntry(tgCfg1[1], "data GeV", "pl");
  leg2->AddEntry(tgCfg1[2], "data Mip (no weights)", "pl");   

  
  TCanvas* chER = new TCanvas();
  chER->cd();
  tgCfg1[0]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg1[0]->GetYaxis()->SetTitle("#sigma(#SigmaE)/ <#SigmaE>");
  tgCfg1[0]->GetYaxis()->SetRangeUser(0., 0.3);
  tgCfg1[0]->Draw("ap");
  tgCfg1[1]->Draw("p, same");
  legTGM->Draw("same");
  chER->Print("dataVsMC/SimGeV_DataGeV_layers8_config1.png", "png");
  chER->Print("dataVsMC/SimGeV_DataGeV_layers8_config1.root", "root");
  //
  TCanvas* chER2 = new TCanvas();
  chER2->cd();
  tgCfg2[0]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg2[0]->GetYaxis()->SetTitle("#sigma(#SigmaE)/ <#SigmaE>");
  tgCfg2[0]->GetYaxis()->SetRangeUser(0., 0.3);
  tgCfg2[0]->Draw("ap");
  tgCfg2[1]->Draw("p, same");
  legTGM2->Draw("same");
  chER2->Print("dataVsMC/SimGeV_DataGeV_layers8_config2.png", "png");
  chER2->Print("dataVsMC/SimGeV_DataGeV_layers8_config2.root", "root");

  //////////
  TCanvas* chEM = new TCanvas();
  chEM->cd();
  tgCfg1[5]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg1[5]->GetYaxis()->SetTitle("<#SigmaE> / beam energy");
  tgCfg1[5]->GetYaxis()->SetRangeUser(0., 1.2);
  tgCfg1[5]->Draw("ap");
  tgCfg1[6]->Draw("p, same");
  legTGM->Draw("same");
  chEM->Print("dataVsMC/scale_SimGeV_DataGeV_layers8_config1.png", "png");
  chEM->Print("dataVsMC/scale_SimGeV_DataGeV_layers8_config1.root", "root");
  //
  TCanvas* chEM2 = new TCanvas();
  chEM2->cd();
  tgCfg2[5]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg2[5]->GetYaxis()->SetTitle("<#SigmaE> beam energy");
  tgCfg2[5]->GetYaxis()->SetRangeUser(0., 1.2);
  tgCfg2[5]->Draw("ap");
  tgCfg2[6]->Draw("p, same");
  legTGM2->Draw("same");
  chEM2->Print("dataVsMC/scale_SimGeV_DataGeV_layers8_config2.png", "png");
  chEM2->Print("dataVsMC/scale_SimGeV_DataGeV_layers8_config2.root", "root");
 
  std::cout << " >>> plot linearity " << std::endl;
  ////////
  TCanvas* chLi = new TCanvas();
  chLi->cd();
  tgCfg1[3]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg1[3]->GetYaxis()->SetTitle("#SigmaE (GeV)");
  //  tgCfg1[3]->GetYaxis()->SetRangeUser(0., 0.3);
  tgCfg1[3]->Draw("ap");
  tgCfg1[4]->Draw("p, same");
  legTGM->Draw("same");
  chLi->Print("dataVsMC/linearity_SimGeV_DataGeV_layers8_config1.png", "png");
  chLi->Print("dataVsMC/linearity_SimGeV_DataGeV_layers8_config1.root", "root");
  //
  TCanvas* chLi2 = new TCanvas();
  chLi2->cd();
  tgCfg2[3]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg2[3]->GetYaxis()->SetTitle("#SigmaE (GeV)");
  //  tgCfg2[3]->GetYaxis()->SetRangeUser(0., 0.3);
  tgCfg2[3]->Draw("ap");
  tgCfg2[4]->Draw("p, same");
  legTGM2->Draw("same");
  chLi2->Print("dataVsMC/linearity_SimGeV_DataGeV_layers8_config2.png", "png");
  chLi2->Print("dataVsMC/linearity_SimGeV_DataGeV_layers8_config2.root", "root");
  ////////



  TCanvas* chDA = new TCanvas();
  chDA->cd();
  tgCfg1[1]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg1[1]->GetYaxis()->SetTitle("#sigma(#SigmaE)/ <#SigmaE>");
  tgCfg1[1]->GetYaxis()->SetRangeUser(0., 0.3);
  tgCfg1[1]->Draw("ap");
  tgCfg1[2]->Draw("p, same");
  leg->Draw("same");
  chDA->Print("dataVsMC/DataMip_DataGeV_layers8_config1.png", "png");
  chDA->Print("dataVsMC/DataMip_DataGeV_layers8_config1.root", "root");

  TCanvas* chDA2 = new TCanvas();
  chDA2->cd();
  tgCfg2[1]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfg2[1]->GetYaxis()->SetTitle("#sigma(#SigmaE)/ <#SigmaE>");
  tgCfg2[1]->GetYaxis()->SetRangeUser(0., 0.3);
  tgCfg2[1]->Draw("ap");
  tgCfg2[2]->Draw("p, same");
  leg2->Draw("same");
  chDA2->Print("dataVsMC/DataMip_DataGeV_layers8_config2.png", "png");
  chDA2->Print("dataVsMC/DataMip_DataGeV_layers8_config2.root", "root");
 

  int iColorsM[5] = {kRed, kCyan, kBlue, kGreen+2, kYellow-1};
  int energyPM[5] = {20, 70, 100, 200, 250};                                                                                                                                   



  TFile* inFM[4];
  inFM[0] = TFile::Open("analyzed_DATA_layers8_config1.root");
  inFM[1] = TFile::Open("analyzed_DATA_layers8_config2.root");
  inFM[2] = TFile::Open("analyzed_SIM_layers8_config1.root");
  inFM[3] = TFile::Open("analyzed_SIM_layers8_config2.root");

  TGraphErrors* tgM_daCfg1[5];
  TGraphErrors* tgM_daCfg2[5];
  TGraphErrors* tgM_mcCfg1[5];
  TGraphErrors* tgM_mcCfg2[5];
  for(int iF=0; iF<5; ++iF){
    tgM_daCfg1[iF] = (TGraphErrors*)inFM[0]->Get(Form("enLayer_MIP_E%d", energyPM[iF]));
    tgM_daCfg2[iF] = (TGraphErrors*)inFM[1]->Get(Form("enLayer_MIP_E%d", energyPM[iF]));
    tgM_mcCfg1[iF] = (TGraphErrors*)inFM[2]->Get(Form("enLayer_MIP_E%d", energyPM[iF]));
    tgM_mcCfg2[iF] = (TGraphErrors*)inFM[3]->Get(Form("enLayer_MIP_E%d", energyPM[iF]));

    tgM_daCfg1[iF]->SetMarkerStyle(20);
    tgM_daCfg2[iF]->SetMarkerStyle(20);
    tgM_mcCfg1[iF]->SetMarkerStyle(21);
    tgM_mcCfg2[iF]->SetMarkerStyle(21);

    tgM_daCfg1[iF]->SetMarkerColor(iColorsM[iF]);
    tgM_daCfg2[iF]->SetMarkerColor(iColorsM[iF]);
    tgM_mcCfg1[iF]->SetMarkerColor(iColorsM[iF]);
    tgM_mcCfg2[iF]->SetMarkerColor(iColorsM[iF]);

  }


  TLegend *legM = new TLegend(0.80,0.70,0.98,0.95,NULL,"brNDC");
  legM->SetTextFont(42);
  legM->SetFillColor(kWhite);
  legM->SetLineColor(kWhite);
  legM->SetShadowColor(kWhite);
  for(int iP=0; iP<5; ++iP){
    legM->AddEntry(tgM_daCfg1[iP], Form("%dGeV", energyPM[iP]), "p");
  }


  TCanvas* cTM = new TCanvas();
  cTM->cd();
  tgM_daCfg1[0]->GetXaxis()->SetTitle("shower depth (X_{0})");
  tgM_daCfg1[0]->GetYaxis()->SetTitle("#Sigma_{layer} E_{layer} (MIP)");
  tgM_daCfg1[0]->GetXaxis()->SetRangeUser(0., 30.);
  tgM_daCfg1[0]->GetYaxis()->SetRangeUser(0., 2.5e3);
  tgM_daCfg1[0]->Draw("ap");
  tgM_mcCfg1[0]->Draw("p, same");
  for(int iT=1; iT<5; ++iT){
    tgM_daCfg1[iT]->Draw("p, same");
    tgM_mcCfg1[iT]->Draw("p, same");
  }
  legM->Draw("same");
  cTM->Print("dataVsMC/enLayer_MIP_nLayer8_config1.png", "png");
  cTM->Print("dataVsMC/enLayer_MIP_nLaye8_config1.root", "root");

  TCanvas* cTM2 = new TCanvas();
  cTM2->cd();
  tgM_daCfg2[0]->GetXaxis()->SetTitle("shower depth (X_{0})");
  tgM_daCfg2[0]->GetYaxis()->SetTitle("#Sigma_{layer} E_{layer} (MIP)");
  tgM_daCfg2[0]->GetXaxis()->SetRangeUser(0., 30.);
  tgM_daCfg2[0]->GetYaxis()->SetRangeUser(0., 2.5e3);
  tgM_daCfg2[0]->Draw("ap");
  tgM_mcCfg2[0]->Draw("p, same");
  for(int iT=1; iT<5; ++iT){
    tgM_daCfg2[iT]->Draw("p, same");
    tgM_mcCfg2[iT]->Draw("p, same");
  }
  legM->Draw("same");
  cTM2->Print("dataVsMC/enLayer_MIP_nLayer8_config2.png", "png");
  cTM2->Print("dataVsMC/enLayer_MIP_nLaye8_config2.root", "root");


}
