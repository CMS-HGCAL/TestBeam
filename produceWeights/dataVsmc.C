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
  int iColors[3] = {kRed, kBlue, kGreen+2};
  //  int energyP[3] = {20, 70, 100, 200, 250};
  
  TFile* inF[4];
  inF[0] = TFile::Open("SIM_resolution_layers8_config1.root");
  inF[1] = TFile::Open("SIM_resolution_layers8_config2.root");
  inF[2] = TFile::Open("DATA_resolution_layers8_config1.root");
  inF[3] = TFile::Open("DATA_resolution_layers8_config2.root");


  TGraphErrors* tgCfg1[3];
  TGraphErrors* tgCfg2[3];

  tgCfg1[0] = (TGraphErrors*)inF[0]->Get("resolution_GeV");
  tgCfg1[1] = (TGraphErrors*)inF[2]->Get("resolution_GeV");
  tgCfg1[2] = (TGraphErrors*)inF[2]->Get("resolution_Mip");

  tgCfg2[0] = (TGraphErrors*)inF[1]->Get("resolution_GeV");
  tgCfg2[1] = (TGraphErrors*)inF[3]->Get("resolution_GeV");
  tgCfg2[2] = (TGraphErrors*)inF[3]->Get("resolution_Mip");


  for(int iT=0; iT<3; ++iT){
    tgCfg1[iT]->SetMarkerStyle(20);
    tgCfg2[iT]->SetMarkerStyle(20);

    tgCfg1[iT]->SetMarkerColor(iColors[iT]);
    tgCfg2[iT]->SetMarkerColor(iColors[iT]);
  }

  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.65,0.65,0.85,0.85,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->SetHeader("6 X_{0} - 15 X_{0}");
  legTGM->AddEntry(tgCfg1[0], "sim GeV", "pl");
  legTGM->AddEntry(tgCfg1[1], "data GeV", "pl");
  //  legTGM->AddEntry(tgCfg1[2], "data Mip", "l");

  TLegend *legTGM2 = new TLegend(0.65,0.65,0.85,0.85,NULL,"brNDC");
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
 
 
}
