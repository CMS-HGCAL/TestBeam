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



void CERNplusFNAL(){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

               
  int nLayers = 8;
  int config = 2;
  int iColors[6] = {kRed, kBlue, kGreen+2, kCyan, kRed, kBlue};
  //  int energyP[3] = {20, 70, 100, 200, 250};
  
  TFile* inF[5];
  inF[0] = TFile::Open("SIM_resolution_layers8_config2.root");
  inF[1] = TFile::Open("DATA_resolution_layers8_config2.root");
  inF[2] = TFile::Open("fromFNAL/Relative_Resolution_BeamEnergy_EMM.root");
  inF[3] = TFile::Open("fromFNAL/Fractional_Response_BeamEnergy_EMM.root");
  inF[4] = TFile::Open("fromFNAL/Response_Linearity_BeamEnergy_EMM.root");


  TGraphErrors* tgCfgR2[6];
  TGraphErrors* tgCfgL2[6];
  TGraphErrors* tgCfgS2[6];
  tgCfgR2[0] = (TGraphErrors*)inF[0]->Get("resolution_GeV");
  tgCfgR2[1] = (TGraphErrors*)inF[1]->Get("resolution_GeV");
  tgCfgR2[2] = (TGraphErrors*)inF[2]->Get("Graph_Sim");
  tgCfgR2[3] = (TGraphErrors*)inF[2]->Get("Graph_Data");

  tgCfgL2[0] = (TGraphErrors*)inF[0]->Get("linerity_GeV");
  tgCfgL2[1] = (TGraphErrors*)inF[1]->Get("linearity_GeV");
  tgCfgL2[2] = (TGraphErrors*)inF[4]->Get("Graph_Sim");
  tgCfgL2[3] = (TGraphErrors*)inF[4]->Get("Graph_Data");

  tgCfgS2[0] = (TGraphErrors*)inF[0]->Get("mean_GeV");
  tgCfgS2[1] = (TGraphErrors*)inF[1]->Get("mean_GeV");
  tgCfgS2[2] = (TGraphErrors*)inF[3]->Get("Graph_Sim");
  tgCfgS2[3] = (TGraphErrors*)inF[3]->Get("Graph_Data");

  tgCfgR2[4] = new TGraphErrors();
  tgCfgR2[5] = new TGraphErrors();
  tgCfgL2[4] = new TGraphErrors();
  tgCfgL2[5] = new TGraphErrors();
  tgCfgS2[4] = new TGraphErrors();
  tgCfgS2[5] = new TGraphErrors();

  for(int iP=0; iP<7; ++iP){
    Double_t xP, yP;
    tgCfgR2[2]->GetPoint(iP, xP, yP);
    tgCfgR2[4]->SetPoint(iP, xP, yP);
    Double_t ey = tgCfgR2[2]->GetErrorY(iP);
    tgCfgR2[4]->SetPointError(iP, 0, ey);

    tgCfgR2[3]->GetPoint(iP, xP, yP);
    tgCfgR2[5]->SetPoint(iP, xP, yP);
    ey = tgCfgR2[3]->GetErrorY(iP);
    tgCfgR2[5]->SetPointError(iP, 0, ey);

    tgCfgL2[2]->GetPoint(iP, xP, yP);
    tgCfgL2[4]->SetPoint(iP, xP, yP/1000.);
    ey = tgCfgL2[2]->GetErrorY(iP);
    tgCfgL2[4]->SetPointError(iP, 0, ey/1000.);
   
    tgCfgL2[3]->GetPoint(iP, xP, yP);
    tgCfgL2[5]->SetPoint(iP, xP, yP/1000.);
    ey = tgCfgL2[3]->GetErrorY(iP);
    tgCfgL2[5]->SetPointError(iP, 0, ey/1000.);

    tgCfgS2[2]->GetPoint(iP, xP, yP);
    tgCfgS2[4]->SetPoint(iP, xP, yP);
    ey = tgCfgS2[2]->GetErrorY(iP);
    tgCfgS2[4]->SetPointError(iP, 0, ey);
   
    tgCfgS2[3]->GetPoint(iP, xP, yP);
    tgCfgS2[5]->SetPoint(iP, xP, yP);
    ey = tgCfgS2[3]->GetErrorY(iP);
    tgCfgS2[5]->SetPointError(iP, 0, ey);
  }
  for(iP=0; iP<4; ++iP){
    Double_t xP, yP;
    tgCfgR2[0]->GetPoint(iP+2, xP, yP);
    tgCfgR2[4]->SetPoint(tgCfgR2[4]->GetN(), xP, yP);
    Double_t ey = tgCfgR2[0]->GetErrorY(iP+2);
    tgCfgR2[4]->SetPointError(tgCfgR2[4]->GetN()-1, 0, ey);
    std::cout << " iP = " << iP << " x = " << xP << " y = " << yP << " eY = " << ey << std::endl;

    tgCfgR2[1]->GetPoint(iP+3, xP, yP);
    tgCfgR2[5]->SetPoint(tgCfgR2[5]->GetN(), xP, yP);
    ey = tgCfgR2[1]->GetErrorY(iP+3);
    tgCfgR2[5]->SetPointError(tgCfgR2[5]->GetN()-1, 0, ey);

    tgCfgL2[0]->GetPoint(iP+2, xP, yP);
    tgCfgL2[4]->SetPoint(tgCfgL2[4]->GetN(), xP, yP);
    ey = tgCfgL2[0]->GetErrorY(iP+2);
    tgCfgL2[4]->SetPointError(tgCfgL2[4]->GetN()-1, 0, ey);
   
    tgCfgL2[1]->GetPoint(iP+3, xP, yP);
    tgCfgL2[5]->SetPoint(tgCfgL2[5]->GetN(), xP, yP);
    ey = tgCfgL2[1]->GetErrorY(iP+3);
    tgCfgL2[5]->SetPointError(tgCfgL2[5]->GetN()-1, 0, ey);

    tgCfgS2[0]->GetPoint(iP+2, xP, yP);
    tgCfgS2[4]->SetPoint(tgCfgS2[4]->GetN(), xP, yP);
    ey = tgCfgS2[0]->GetErrorY(iP+2);
    tgCfgS2[4]->SetPointError(tgCfgS2[4]->GetN()-1, 0, ey);
   
    tgCfgS2[1]->GetPoint(iP+3, xP, yP);
    tgCfgS2[5]->SetPoint(tgCfgS2[5]->GetN(), xP, yP);
    ey = tgCfgS2[1]->GetErrorY(iP+3);
    tgCfgS2[5]->SetPointError(tgCfgS2[5]->GetN()-1, 0, ey);
    //    std::cout << "scale  iP = " << iP << " x = " << xP << " y = " << yP << " eY = " << ey << std::endl;
    ey = tgCfgS2[5]->GetErrorY(tgCfgS2[5]->GetN()-1);
    //    std::cout << "scale  iP = " << iP << " x = " << xP << " y = " << yP << " eY = " << ey << std::endl;
  }
  // tgCfgR2[4]->SetPoint(tgCfgR2[4]->GetN(), 400, 5);
  // tgCfgR2[5]->SetPoint(tgCfgR2[4]->GetN(), 400, 5);
  // tgCfgS2[4]->SetPoint(tgCfgS2[4]->GetN(), 400, 5);
  // tgCfgS2[5]->SetPoint(tgCfgS2[4]->GetN(), 400, 5);
  // tgCfgL2[4]->SetPoint(tgCfgR2[4]->GetN(), 400, 400);
  // tgCfgL2[5]->SetPoint(tgCfgR2[4]->GetN(), 400, 400);
  
  for(int iT=0; iT<6; ++iT){
    tgCfgR2[iT]->SetMarkerStyle(20);
    tgCfgR2[iT]->SetMarkerSize(1);
    tgCfgR2[iT]->SetMarkerColor(iColors[iT]);
    tgCfgR2[iT]->SetLineColor(iColors[iT]);
    tgCfgR2[iT]->SetLineWidth(1);

    tgCfgL2[iT]->SetMarkerStyle(20);
    tgCfgL2[iT]->SetMarkerSize(1);
    tgCfgL2[iT]->SetMarkerColor(iColors[iT]);
    tgCfgL2[iT]->SetLineColor(iColors[iT]);
    tgCfgL2[iT]->SetLineWidth(1);

    tgCfgS2[iT]->SetMarkerStyle(20);
    tgCfgS2[iT]->SetMarkerSize(1);
    tgCfgS2[iT]->SetMarkerColor(iColors[iT]);
    tgCfgS2[iT]->SetLineColor(iColors[iT]);
    tgCfgS2[iT]->SetLineWidth(1);
  }

  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.65,0.65,0.85,0.85,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  // legTGM->AddEntry(tgCfgR2[0], "CERN sim", "pl");
  // legTGM->AddEntry(tgCfgR2[1], "CERN data", "pl");
  // legTGM->AddEntry(tgCfgR2[2], "FNAL sim", "pl");
  // legTGM->AddEntry(tgCfgR2[3], "FNAL data", "pl");
  legTGM->AddEntry(tgCfgR2[4], "sim combined", "pl");
  legTGM->AddEntry(tgCfgR2[5], "data combined", "pl");
  //  legTGM->AddEntry(tgCfg1[2], "data Mip", "l");

  
  gStyle->SetOptFit(1);
  TF1* fitReso = new TF1("fitReso", "sqrt( pow([0]/sqrt(x), 2.) + pow([1]/x, 2.) + pow([2], 2.))", 2., 300.);
  fitReso->SetParName(0, "S");
  fitReso->SetParName(1, "N");
  fitReso->SetParName(2, "C");
  fitReso->SetParLimits(1, 0., 0.16);

  //  tgCfgR2[4]->Fit("fitReso", "R");
  fitReso->SetLineColor(kRed);
  fitReso->SetLineWidth(2);

  TF1* fitReso_D = new TF1("fitReso_D", "sqrt( pow([0]/sqrt(x), 2.) + pow([1]/x, 2.) + pow([2], 2.))", 2., 300.);
  fitReso_D->SetParName(0, "S");
  fitReso_D->SetParName(1, "N");
  fitReso_D->SetParName(2, "C");
  fitReso_D->SetParLimits(1, 0., 0.16);

  //  tgCfgR2[5]->Fit("fitReso_D", "R");
  fitReso_D->SetLineColor(kBlue);
  fitReso_D->SetLineWidth(2);

  TCanvas* chER = new TCanvas();
  chER->cd();
  tgCfgR2[4]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfgR2[4]->GetYaxis()->SetTitle("#sigma(#SigmaE)/ <#SigmaE>");
  tgCfgR2[4]->GetXaxis()->SetRangeUser(0., 300.);
  tgCfgR2[4]->GetYaxis()->SetRangeUser(0., 0.3);
  tgCfgR2[4]->Draw("ap");
  tgCfgR2[5]->Draw("p, same");
  // tgCfgR2[0]->Draw("p, same");
  // tgCfgR2[1]->Draw("p, same");
  // tgCfgR2[2]->Draw("p, same");
  // tgCfgR2[3]->Draw("p, same");
  //  fitReso_D->Draw("same");
  //  fitReso->Draw("same");
  legTGM->Draw("same");
  chER->Print("CERNFNAL/resolution.png", "png");
  chER->Print("CERNFNAL/resolution.root", "root");

  //
  TCanvas* chES = new TCanvas();
  chES->cd();
  tgCfgS2[4]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfgS2[4]->GetYaxis()->SetTitle("<#SigmaE>/ beam energy");
  tgCfgS2[4]->GetXaxis()->SetRangeUser(0., 300.);
  tgCfgS2[4]->GetYaxis()->SetRangeUser(0., 1.1);
  tgCfgS2[4]->Draw("ap");
  tgCfgS2[5]->Draw("p, same");
  // tgCfgS2[0]->Draw("p, same");
  // tgCfgS2[1]->Draw("p, same");
  // tgCfgS2[2]->Draw("p, same");
  // tgCfgS2[3]->Draw("p, same");
  legTGM->Draw("same");
  chES->Print("CERNFNAL/scale.png", "png");
  chES->Print("CERNFNAL/scale.root", "root");


  TCanvas* chEL = new TCanvas();
  chEL->cd();
  tgCfgL2[4]->GetXaxis()->SetTitle("beam energy GeV");
  tgCfgL2[4]->GetYaxis()->SetTitle("<#SigmaE> GeV");
  tgCfgL2[4]->GetXaxis()->SetRangeUser(0., 300.);
  tgCfgL2[4]->GetYaxis()->SetRangeUser(0., 300.);
  tgCfgL2[4]->Draw("ap");
  tgCfgL2[5]->Draw("p, same");
  // tgCfgL2[0]->Draw("p, same");
  // tgCfgL2[1]->Draw("p, same");
  // tgCfgL2[2]->Draw("p, same");
  // tgCfgL2[3]->Draw("p, same");
  legTGM->Draw("same");
  chEL->Print("CERNFNAL/linearity.png", "png");
  chEL->Print("CERNFNAL/linearity.root", "root");


}
