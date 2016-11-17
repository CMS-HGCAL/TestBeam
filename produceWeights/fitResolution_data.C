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



void fitResolution_data(){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  
  /*
  int nLayers = 28;
  int config = 1;
  int iColors[3] = {kRed, kBlue, kGreen+2};
  int energyP[3] = {20, 100, 200};
  */


  
  int nLayers = 8;
  int config = 1;
  int iColors[6] = {kRed, kViolet, kCyan, kBlue, kGreen+2, kYellow-1};
  int energyP[6] = {20, 32, 70, 100, 200, 250};
  

  /*                          
  int nLayers = 8;
  int config = 2;
  int iColors[6] = {kRed, kViolet, kCyan, kBlue, kGreen+2, kYellow-1}; //, kBlack};
  int energyP[6] = {20, 32, 70, 100, 200, 250};
  */


  TGraphErrors* tg[6];
  TGraphErrors* tgL[6];
  TGraphErrors* tgR[6];
  TGraphErrors* tgS[6];
  for(int iT=0; iT<6; ++iT){
    if(nLayers == 28 && iT >= 3) continue;
    tg[iT] = new TGraphErrors();
    tg[iT]->SetName(Form("mean_E%d", energyP[iT]));
    tg[iT]->SetMarkerStyle(20);

    tgL[iT] = new TGraphErrors();
    tgL[iT]->SetName(Form("linearity_E%d", energyP[iT]));
    tgL[iT]->SetMarkerStyle(20);

    tgR[iT] = new TGraphErrors();
    tgR[iT]->SetName(Form("resolution_Mip%d", energyP[iT]));
    tgR[iT]->SetMarkerStyle(20);

    tgS[iT] = new TGraphErrors();
    tgS[iT]->SetName(Form("resolutionS_E%d", energyP[iT]));
    tgS[iT]->SetMarkerStyle(20);
  }

  TFile* inF[5];
  if(nLayers == 8 && config == 1){
    for(int iP=0; iP<6; ++iP)
      inF[iP] = TFile::Open(Form("showers_ROOT_data/analyzed_Ele%d.root", energyP[iP]));
  }
  else if(nLayers == 8 && config == 2){
    for(int iP=0; iP<6; ++iP)
      inF[iP] = TFile::Open(Form("showers_ROOT_data/analyzed_endShower_Ele%d.root", energyP[iP]));
  }
  else{
    for(int iP=0; iP<3; ++iP)
      inF[iP] = TFile::Open(Form("showers_ROOT_data/analyzed_28layers_Ele%d.root", energyP[iP]));
  }


   TH1F* recoEnergy[6];
   TH1F* recoEnergy_up[6];
   TH1F* recoEnergy_dw[6];

   TH1F* recoMip[6];
   TH1F* recoMip_up[6];
   TH1F* recoMip_dw[6];

  for(int iC=0; iC<6; ++iC){
    if(nLayers == 28 && iC >=3) continue;
    recoEnergy[iC] = (TH1F*)(inF[iC]->Get("h_eAll_all_AbsW_GeV"));
    recoEnergy_up[iC] = (TH1F*)(inF[iC]->Get("h_eAll_all_AbsW_GeV_up"));
    recoEnergy_dw[iC] = (TH1F*)(inF[iC]->Get("h_eAll_all_AbsW_GeV_dw"));
    
    recoEnergy[iC]->SetLineColor(iColors[iC]);
    recoEnergy[iC]->SetLineWidth(2);

    recoMip[iC] = (TH1F*)(inF[iC]->Get("h_eAll_all"));
    recoMip_up[iC] = (TH1F*)(inF[iC]->Get("h_eAll_all_up"));
    recoMip_dw[iC] = (TH1F*)(inF[iC]->Get("h_eAll_all_dw"));

    recoMip[iC]->SetLineColor(iColors[iC]);
    recoMip[iC]->SetLineWidth(2);

    
    recoEnergy[iC]->Rebin(2);
    recoEnergy_up[iC]->Rebin(2);
    recoEnergy_dw[iC]->Rebin(2);

    recoMip[iC]->Rebin(10);
    recoMip_up[iC]->Rebin(10);
    recoMip_dw[iC]->Rebin(10);
    
  }



  TF1* fitF[6];
  TF1* fitF_up[6];
  TF1* fitF_dw[6];

  TF1* fitM[6];
  TF1* fitM_up[6];
  TF1* fitM_dw[6];
  for(int iP=0; iP<6; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    fitF[iP] = new TF1(Form("E%d",energyP[iP]), "gaus", recoEnergy[iP]->GetMean() - 0.5* recoEnergy[iP]->GetRMS(), recoEnergy[iP]->GetMean() + 2.* recoEnergy[iP]->GetRMS());
    fitF_up[iP] = new TF1(Form("E%d_up",energyP[iP]), "gaus", recoEnergy_up[iP]->GetMean() - 0.5* recoEnergy_up[iP]->GetRMS(), recoEnergy_up[iP]->GetMean() + 2.* recoEnergy_up[iP]->GetRMS());
    fitF_dw[iP] = new TF1(Form("E%d_dw",energyP[iP]), "gaus", recoEnergy_dw[iP]->GetMean() - 0.5* recoEnergy_dw[iP]->GetRMS(), recoEnergy_dw[iP]->GetMean() + 2.* recoEnergy_dw[iP]->GetRMS());
    fitM[iP] = new TF1(Form("M%d",energyP[iP]), "gaus", recoMip[iP]->GetMean() - 0.5* recoMip[iP]->GetRMS(), recoMip[iP]->GetMean() + 2.* recoMip[iP]->GetRMS());
    fitM_up[iP] = new TF1(Form("M%d_up",energyP[iP]), "gaus", recoMip_up[iP]->GetMean() - 0.5* recoMip_up[iP]->GetRMS(), recoMip_up[iP]->GetMean() + 2.* recoMip_up[iP]->GetRMS());
    fitM_dw[iP] = new TF1(Form("M%d_dw",energyP[iP]), "gaus", recoMip_dw[iP]->GetMean() - 0.5* recoMip_dw[iP]->GetRMS(), recoMip_dw[iP]->GetMean() + 2.* recoMip_dw[iP]->GetRMS());

    fitF[iP]->SetParameters(recoEnergy[iP]->Integral(), recoEnergy[iP]->GetMean(), recoEnergy[iP]->GetRMS());
    fitF[iP]->SetLineColor(iColors[iP]);

    fitM[iP]->SetParameters(recoMip[iP]->Integral(), recoMip[iP]->GetMean(), recoMip[iP]->GetRMS());
    fitM[iP]->SetLineColor(iColors[iP]);
    //    std::out << << std::endl;
  }

  std::cout << " ci sono " << std::endl;
  //  return;
  TLegend *legTGM = new TLegend(0.80,0.30,0.95,0.55,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iP=0; iP<6; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    legTGM->AddEntry(recoEnergy[iP], Form("%dGeV", energyP[iP]), "l");
  }


  TLegend *leg = new TLegend(0.80,0.70,0.98,0.95,NULL,"brNDC");
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  for(int iP=0; iP<6; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    leg->AddEntry(recoEnergy[iP], Form("%dGeV", energyP[iP]), "l");
  }


   std::string folder = "plots_Energy_data";

   /////////////////////////////

  for(int iP=0; iP<6; ++iP){
    if(nLayers == 28 && iP >=3) continue;
    recoEnergy[iP]->Fit(Form("E%d", energyP[iP]), "RQ");
    recoEnergy_up[iP]->Fit(Form("E%d_up", energyP[iP]), "RQ");
    recoEnergy_dw[iP]->Fit(Form("E%d_dw", energyP[iP]), "RQ");

    recoMip[iP]->Fit(Form("M%d", energyP[iP]), "RQ");
    recoMip_up[iP]->Fit(Form("M%d_up", energyP[iP]), "RQ");
    recoMip_dw[iP]->Fit(Form("M%d_dw", energyP[iP]), "RQ");

    //    std::cout << "*********** mean fit = " << fitM[iP]->GetParameter(1) << std::endl;
  }

  // return;
  
   TCanvas* chE = new TCanvas();
   chE->cd();
   recoEnergy[0]->GetXaxis()->SetTitle("#Sigma E_{i} (GeV)");
   recoEnergy[0]->GetXaxis()->SetRangeUser(0.7, 1.3);
   recoEnergy[0]->GetYaxis()->SetRangeUser(0., 100.e+03);
   if(nLayers == 8){
     recoEnergy[0]->GetXaxis()->SetRangeUser(0., 300.);
     recoEnergy[0]->GetYaxis()->SetRangeUser(0., 1.2e+03);
   }
   recoEnergy[0]->Draw();
   for(int iP=1; iP<6; ++iP){
     if(nLayers == 28 && iP >=3) continue;
     recoEnergy[iP]->Draw("same");
   }
   leg->Draw("same");
   chE->Print(Form((folder+"/energyPlots_layers%d_config%d.png").c_str(), nLayers, config), "png");
   chE->Print(Form((folder+"/energyPlots_layers%d_config%d.root").c_str(), nLayers, config), "root");

   TCanvas* chM = new TCanvas();
   chM->cd();
   recoMip[0]->GetXaxis()->SetTitle("#Sigma E_{i} (MIP)");
   recoMip[0]->GetXaxis()->SetRangeUser(0.7, 1.3);
   recoMip[0]->GetYaxis()->SetRangeUser(0., 100.e+03);
   if(nLayers == 8){
     recoMip[0]->GetXaxis()->SetRangeUser(0., 15.e+03);
     recoMip[0]->GetYaxis()->SetRangeUser(0., 1.e+03);
   }
   recoMip[0]->Draw();
   for(int iP=1; iP<6; ++iP){
     if(nLayers == 28 && iP >=3) continue;
     recoMip[iP]->Draw("same");
   }
   leg->Draw("same");
   chM->Print(Form((folder+"/mipPlots_layers%d_config%d.png").c_str(), nLayers, config), "png");
   chM->Print(Form((folder+"/mipPlots_layers%d_config%d.root").c_str(), nLayers, config), "root");

   //   return;


   for(int iT=0; iT<6; ++iT){
     if(nLayers == 28 && iT >=3) continue;

     tg[0]->SetPoint(iT+1, energyP[iT], fitF[iT]->GetParameter(1)/energyP[iT]);
     tg[0]->SetPointError(iT+1, 0., sqrt(pow(fitF[iT]->GetParError(1)/energyP[iT], 2.) + pow(abs(fitF_up[iT]->GetParameter(1) - fitF_dw[iT]->GetParameter(1)) /2. / energyP[iT], 2.) ) );


     tgL[0]->SetPoint(iT+1, energyP[iT], fitF[iT]->GetParameter(1));
     tgL[0]->SetPointError(iT+1, 0., sqrt(pow(fitF[iT]->GetParError(1), 2.) + pow(abs(fitF_up[iT]->GetParameter(1) - fitF_dw[iT]->GetParameter(1)) /2., 2.) ) );

     std::cout << " >>>>>>> in MIP energy = " << energyP[iT] << " reso = " << fitM[iT]->GetParameter(2) / fitM[iT]->GetParameter(1)  << std::endl;
     tgR[0]->SetPoint(iT+1, energyP[iT], fitM[iT]->GetParameter(2) / fitM[iT]->GetParameter(1));
     tgR[0]->SetPointError(iT+1, 0., sqrt( pow(fitM[iT]->GetParameter(2) / fitM[iT]->GetParameter(1) * sqrt( pow(fitM[iT]->GetParError(2)/fitM[iT]->GetParameter(2), 2.) + 
       													     pow(fitM[iT]->GetParError(1)/fitM[iT]->GetParameter(1), 2.) ), 2.) +
					   pow(std::abs(fitM_up[iT]->GetParameter(2) / fitM[iT]->GetParameter(1) - fitM_dw[iT]->GetParameter(2) / fitM[iT]->GetParameter(1))/2., 2.) ) );
     //     tgR[0]->SetPointError(iT+1, 0., 100.*pow(fitF[iT]->GetParError(2),2.)/ (1.* energyP[iT]) );
     
     std::cout << " >>>>>>> in GeV energy = " << energyP[iT] << " reso = " << fitF[iT]->GetParameter(2) / fitF[iT]->GetParameter(1)  << std::endl;
     tgS[0]->SetPoint(iT+1, energyP[iT], fitF[iT]->GetParameter(2) / fitF[iT]->GetParameter(1));
     tgS[0]->SetPointError(iT+1, 0., sqrt( pow(fitF[iT]->GetParameter(2) / fitF[iT]->GetParameter(1) * sqrt( pow(fitF[iT]->GetParError(2)/fitF[iT]->GetParameter(2), 2.) + 
      													     pow(fitF[iT]->GetParError(1)/fitF[iT]->GetParameter(1), 2.) ), 2.) +
					   pow(std::abs(fitF_up[iT]->GetParameter(2) / fitF[iT]->GetParameter(1) - fitF_dw[iT]->GetParameter(2) / fitF[iT]->GetParameter(1))/2., 2.) ) );
     //     tgS[0]->SetPointError(iT+1, 0.,sqrt( pow(fitF[iT]->GetParError(2), 2.) + pow(fitF[iT]->GetParError(1), 2.) ) );
     // tgS[0]->SetPointError(iT+1, 0., sqrt( pow(fitF[iT]->GetParError(2)/fitF[iT]->GetParameter(1), 2) + 
     //  					   pow(fitF[iT]->GetParameter(2)/fitF[iT]->GetParameter(1)/fitF[iT]->GetParameter(1) * fitF[iT]->GetParError(1), 2) +
     // 					   pow(std::abs(fitF_up[iT]->GetParameter(2) / fitF_up[iT]->GetParameter(1) - fitF_dw[iT]->GetParameter(2) / fitF_dw[iT]->GetParameter(1)) /2., 2)) );


     //std::cout << " >>> energyP[iT] = " << energyP[iT] << " fitF[iT]->GetParameter(1) / (1.* energyP[iT]) = " << fitF[iT]->GetParameter(1) / (1.* energyP[iT]) << std::endl;
     //     std::cout << " >>> energyP[iT] = " << energyP[iT] << " Y = " << pow(fitF[iT]->GetParameter(2),2) / (1.* energyP[iT]) * 100. << std::endl;

     tg[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgL[0]->SetPoint(tg[iT]->GetN(), 400, 300);
     tgR[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgS[0]->SetPoint(tg[iT]->GetN(), 400, 5);
  }



    /////////TGraph linearity
    TCanvas* tgLi = new TCanvas();
    tgLi->cd();
    tgL[0]->GetYaxis()->SetTitle("#Sigma E_{i} (GeV)");
    tgL[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tgL[0]->GetXaxis()->SetRangeUser(0, 300.);
    if(nLayers == 8) tg[0]->GetYaxis()->SetRangeUser(0., 1.01);
    tgL[0]->Draw("ap");
    tgLi->Print(Form((folder+"/energy_LinearityVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");
    tgLi->Print(Form((folder+"/energy_LinearityVsEnergy_layers%d_config%d.root").c_str(), nLayers, config), "root");

    /////////TGraph
    TCanvas* tgM = new TCanvas();
    tgM->cd();
    tg[0]->GetYaxis()->SetTitle("#Sigma E_{i} / E_{beam}");
    tg[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tg[0]->GetXaxis()->SetRangeUser(0, 300.);
    tg[0]->GetYaxis()->SetRangeUser(0.9, 1.01);
    if(nLayers == 8) tg[0]->GetYaxis()->SetRangeUser(0., 1.01);
    tg[0]->Draw("ap");
    tgM->Print(Form((folder+"/energy_MeanVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");
    tgM->Print(Form((folder+"/energy_MeanVsEnergy_layers%d_config%d.root").c_str(), nLayers, config), "root");


    TFile outData(Form("DATA_resolution_layers%d_config%d.root", nLayers, config), "recreate");
    outData.cd();
    tg[0]->Write("mean_GeV");
    tgL[0]->Write("linearity_GeV");
    tgS[0]->Write("resolution_GeV");
    tgR[0]->Write("resolution_Mip");
    outData.Close();



    gStyle->SetOptFit(1);
    /////////TGraph

    
    TF1* fitResoM;

    //fit 3par all
    //    int type = 1;

    //fit 3par tail
    //int type = 2;

    //fit 2par all
    int type = 3;

    if(type == 1){
      fitResoM = new TF1("fitResoM", "sqrt( pow([0]/sqrt(x), 2.) + pow([1]/x, 2.) + pow([2], 2.))", 10., 300.);
      fitResoM->SetParName(0, "S");
      fitResoM->SetParName(1, "N");
      fitResoM->SetParName(2, "C");
      fitResoM->SetParLimits(1, 0., 0.16);
      //    fitResoM->SetParLimits(2, 0.0001, 0.05);
    }
    
    if(type == 3){
      fitResoM = new TF1("fitResoM", "sqrt( pow([0]/sqrt(x), 2.) + pow([1], 2.))", 10., 300.);
      fitResoM->SetParName(0, "S");
      fitResoM->SetParName(1, "C");
      //    fitResoM->SetParLimits(1, 0., 0.05);
      //    fitResoM->SetParLimits(1, 0.0001, 0.05);
    }

    if(type == 2){
      fitResoM = new TF1("fitResoM", "sqrt( pow([0]/sqrt(x), 2.) + pow([1]/x, 2.) + pow([2], 2.))", 50., 300.);
      fitResoM->SetParName(0, "S");
      fitResoM->SetParName(1, "N");
      fitResoM->SetParName(2, "C");
      fitResoM->SetParLimits(1, 0., 0.16);      
    }

    tgR[0]->Fit("fitResoM", "R");
    fitResoM->SetLineColor(kRed);
    fitResoM->SetLineWidth(2);

    /*
    for(int point=0; point<tgR[0]->GetN(); ++point){
      double ey = tgR[0]->GetErrorY(point);
      tgR[0]->SetPointError(point, 0., ey*sqrt(fitResoM->GetChisquare()/fitResoM->GetNDF()));
      }
    */




    TCanvas* tgResoM = new TCanvas();
    tgResoM->cd();
    tgR[0]->GetYaxis()->SetTitle("#sigma( #Sigma E) / E_{reco}");
    tgR[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tgR[0]->GetXaxis()->SetRangeUser(0, 300.);
    tgR[0]->GetYaxis()->SetRangeUser(0., 0.3);
    tgR[0]->Draw("ap");
    fitResoM->Draw("same");
    //fitResoM_post->Draw("same");
    tgResoM->Print(Form((folder+"/energyMIPnonW_ResolutionVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");
    tgResoM->Print(Form((folder+"/energyMIPnonW_ResolutionVsEnergy_layers%d_config%d.root").c_str(), nLayers, config), "root");



    /////////////////////////////////////////////////////////
    TF1* fitReso;

    //fit 3par all
    //    type = 1;

    //fit 3par tail  
    //    int type = 2;
    //fit 2par all 
    //    int type = 3; 

    if(type == 1){
      fitReso = new TF1("fitReso", "sqrt( pow([0]/sqrt(x), 2.) + + pow([1]/x, 2.) + pow([2], 2.))", 10., 300.);
      fitReso->SetParName(0, "S");
      fitReso->SetParName(1, "N");
      fitReso->SetParName(2, "C");
      fitReso->SetParLimits(1, 0., 0.16);
      //    fitReso->SetParLimits(2, 0.0001, 0.05);
    }

    if(type == 3){
      fitReso = new TF1("fitReso", "sqrt( pow([0]/sqrt(x), 2.) + pow([1], 2.))", 10., 300.);
      fitReso->SetParName(0, "S");
      fitReso->SetParName(1, "C");
      //    fitReso->SetParLimits(1, 0., 0.05);
      fitReso->SetParLimits(1, 0.0001, 0.05);
    }

    if(type == 2){
      fitReso = new TF1("fitReso", "sqrt( pow([0]/sqrt(x), 2.) + pow([1]/x, 2.) + pow([2], 2.))", 50., 300.);
      fitReso->SetParName(0, "S");
      fitReso->SetParName(1, "N");
      fitReso->SetParName(2, "C");
      fitReso->SetParLimits(1, 0., 0.16);
    }
    
    tgS[0]->Fit("fitReso", "R");
    fitReso->SetLineColor(kGreen+2);
    fitReso->SetLineWidth(2);
    
    /*
    for(int point=0; point<tgS[0]->GetN(); ++point){
      double ey = tgS[0]->GetErrorY(point);
      tgS[0]->SetPointError(point,0., ey*sqrt(fitReso->GetChisquare()/fitReso->GetNDF()));
	//	tgS[0]->SetPointEYlow(point,ey*sqrt(fitReso->GetChisquare()/fitReso->GetNDF()));
    }
    */


    /////////TGraph   simple
    TCanvas* tgResoS = new TCanvas();
    tgResoS->cd();
    tgS[0]->GetYaxis()->SetTitle("#sigma( #Sigma E )/ E_{reco}");
    tgS[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tgS[0]->GetXaxis()->SetRangeUser(0., 300.);
    tgS[0]->GetYaxis()->SetRangeUser(0., 0.3);
    tgS[0]->Draw("ap");
    fitReso->Draw("same");
    //    fitReso_post->Draw("same");
    tgResoS->Print(Form((folder+"/energy_SResolutionVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");
    tgResoS->Print(Form((folder+"/energy_SResolutionVsEnergy_layers%d_config%d.root").c_str(), nLayers, config), "root");






    //        return;
    ////////////////
    std::cout << " >>> fatto " << std::endl;
    //    return;
    TGraphErrors* tM[6];
    for(int iT=0; iT<6; ++iT){
      if(nLayers == 28 && iT >= 3) continue;
      tM[iT] = (TGraphErrors*)inF[iT]->Get("enLayer_MIP");
      tM[iT]->SetName(Form("enLayer_MIP_E%d", energyP[iT]));
      
      tM[iT]->SetMarkerStyle(20);
      tM[iT]->SetMarkerColor(iColors[iT]);
      tM[iT]->SetLineColor(iColors[iT]);
      tM[iT]->SetLineWidth(2);
      tM[iT]->SetLineStyle(iT);
    }
    
    TFile newOUT(Form("analyzed_DATA_layers%d_config%d.root", nLayers, config), "recreate");
    newOUT.cd();
    for(int iT=0; iT<6; ++iT){
      tM[iT]->SetName(Form("enLayer_MIP_E%d", energyP[iT]));
      tM[iT]->Write(Form("enLayer_MIP_E%d", energyP[iT]));
    }
    newOUT.Close();


    TCanvas* cTM = new TCanvas();
    cTM->cd();
    tM[0]->GetXaxis()->SetTitle("shower depth (X_{0})");
    tM[0]->GetYaxis()->SetTitle("#Sigma_{layer} E_{layer} (MIP)");
    tM[0]->GetXaxis()->SetRangeUser(0., 300.);
    tM[0]->GetYaxis()->SetRangeUser(0., 2.e3);
    tM[0]->Draw("ap");
    for(int iT=1; iT<6; ++iT){
      if(nLayers == 28 && iT >= 3) continue;
      tM[iT]->Draw("p, same");
    }
    leg->Draw("same");
    cTM->Print(Form((folder+"/enLayer_MIP_nLayer%d_config%d.png").c_str(), nLayers, config), "png");
    cTM->Print(Form((folder+"/enLayer_MIP_nLayer%d_config%d.root").c_str(), nLayers, config), "root");
    

}
