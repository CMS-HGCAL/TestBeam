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



void fitResolution(){
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  
  /*      
  int nLayers = 28;
  int config = 1;
  int iColors[5] = {kRed, kCyan, kBlue, kGreen+2, kYellow-1};
  int energyP[5] = {20, 70, 100, 200, 250};
  */  


  int nLayers = 8;
  int config = 1;
  int iColors[5] = {kRed, kCyan, kBlue, kGreen+2, kYellow-1};
  int energyP[5] = {20, 70, 100, 200, 250};

          
  /*      
  int nLayers = 8;
  int config = 2;
  int iColors[5] = {kRed, kCyan, kBlue, kGreen+2, kYellow-1};
  int energyP[5] = {20, 70, 100, 200, 250};
  */

  TGraphErrors* tg[5];
  TGraphErrors* tgL[5];
  TGraphErrors* tgR[5];
  TGraphErrors* tgS[5];
  for(int iT=0; iT<5; ++iT){
    //    if(nLayers == 28 && iT >= 4) continue;
    tg[iT] = new TGraphErrors();
    tg[iT]->SetName(Form("mean_E%d", energyP[iT]));
    tg[iT]->SetMarkerStyle(20);

    tgL[iT] = new TGraphErrors();
    tgL[iT]->SetName(Form("linearity_E%d", energyP[iT]));
    tgL[iT]->SetMarkerStyle(20);

    tgR[iT] = new TGraphErrors();
    tgR[iT]->SetName(Form("resolution_E%d", energyP[iT]));
    //    tgR[iT]->SetPoint(0, -1, -1);

    //    tgR[iT]->SetLineColor(iColors[iT]);
    //    tgR[iT]->SetMarkerColor(iColors[iT]);
    tgR[iT]->SetMarkerStyle(20);

    tgS[iT] = new TGraphErrors();
    tgS[iT]->SetName(Form("resolutionS_E%d", energyP[iT]));
    //    tgS[iT]->SetPoint(0, -1, -1);

    //    tgS[iT]->SetLineColor(iColors[iT]);
    //    tgS[iT]->SetMarkerColor(iColors[iT]);
    tgS[iT]->SetMarkerStyle(20);
  }


  TFile* inF[5];
  if(nLayers == 8 && config == 1){
    for(int iP=0; iP<5; ++iP)
      inF[iP] = TFile::Open(Form("energyLayer_ROOT/outF_Ele%d.root", energyP[iP]));
  }
  else if(nLayers == 8 && config == 2){
    for(int iP=0; iP<5; ++iP)
      inF[iP] = TFile::Open(Form("energyLayer_ROOT/outF_endShower_Ele%d.root", energyP[iP]));
  }
  else{
    for(int iP=0; iP<5; ++iP)
      inF[iP] = TFile::Open(Form("energyLayer_ROOT/outF_28layers_Ele%d.root", energyP[iP]));
  }

   TH1F* recoEnergy[5];
   TH1F* recoEnergyRel[5];


  for(int iC=0; iC<5; ++iC){
    //    if(nLayers == 28 && iC >=4) continue;
    recoEnergy[iC] = (TH1F*)(inF[iC]->Get("recoEnergy"));
    recoEnergyRel[iC] = (TH1F*)(inF[iC]->Get("recoEnergyRel"));
    
    recoEnergy[iC]->SetLineColor(iColors[iC]);
    recoEnergy[iC]->SetLineWidth(2);

    recoEnergyRel[iC]->SetLineColor(iColors[iC]);
    recoEnergyRel[iC]->SetLineWidth(2);
  }



  TF1* fitF[5];
  TF1* fitFRel[5];
  for(int iP=0; iP<5; ++iP){
    //    if(nLayers == 28 && iP >=4) continue;
    fitF[iP] = new TF1(Form("E%d",energyP[iP]), "gaus", recoEnergy[iP]->GetMean() - 0.5 * recoEnergy[iP]->GetRMS(), recoEnergy[iP]->GetMean() + 2.*recoEnergy[iP]->GetRMS());
    fitFRel[iP] = new TF1(Form("E%dR",energyP[iP]), "gaus", recoEnergy[iP]->GetMean() - 0.5 * recoEnergy[iP]->GetRMS(), recoEnergy[iP]->GetMean() + 2.*recoEnergy[iP]->GetRMS());

    fitF[iP]->SetParameters(recoEnergy[iP]->Integral(), recoEnergy[iP]->GetMean(), recoEnergy[iP]->GetRMS());
    fitF[iP]->SetLineColor(iColors[iP]);

    fitFRel[iP]->SetParameters(recoEnergyRel[iP]->Integral(), recoEnergy[iP]->GetMean(), recoEnergy[iP]->GetRMS());
    fitFRel[iP]->SetLineColor(iColors[iP]);
  }

  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.80,0.30,0.95,0.55,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iP=0; iP<5; ++iP){
    //    if(nLayers == 28 && iP >=4) continue;
    legTGM->AddEntry(recoEnergy[iP], Form("%dGeV", energyP[iP]), "l");
  }


  TLegend *leg = new TLegend(0.80,0.70,0.98,0.95,NULL,"brNDC");
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  for(int iP=0; iP<5; ++iP){
    //    if(nLayers == 28 && iP >=4) continue;
    leg->AddEntry(recoEnergy[iP], Form("%dGeV", energyP[iP]), "l");
  }


  TLatex t1;
  t1.SetNDC();
  t1.SetTextSize(0.03);
  t1.SetTextFont(132);
  t1.SetTextColor(kRed);

  TLatex t2;
  t2.SetNDC();
  t2.SetTextSize(0.03);
  t2.SetTextFont(132);
  t2.SetTextColor(kBlue);


  TLatex t3;
  t3.SetNDC();
  t3.SetTextSize(0.03);
  t3.SetTextFont(132);
  t3.SetTextColor(kGreen+2);

  for(int iP=0; iP<5; ++iP){
    //    if(nLayers == 28 && iP >=4) continue;
    recoEnergy[iP]->Fit(Form("E%d", energyP[iP]), "R");
  }

  //    float energyValues[5] = {10., 30., 60., 100., 250.};

  //  std::string folder = "plots"+std::string(Form("_E%d",ene));
  std::string folder = "plots_Energy";
  //     std::string folder = "plots_E250";
  
  /*
  TCanvas* chE = new TCanvas();
  chE->cd();
  recoEnergy[0]->GetXaxis()->SetTitle("#Sigma E_{i}");
  recoEnergy[0]->Draw();
  for(int iP=1; iP<5; ++iP){
    recoEnergy[iP]->Draw("same");
  }
   t1.DrawLatex(0.2,0.9,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF[0]->GetParameter(1), fitF[0]->GetParError(1), 
   			     fitF[0]->GetParameter(2), fitF[0]->GetParError(2)));
   t2.DrawLatex(0.2,0.85,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF[1]->GetParameter(1), fitF[1]->GetParError(1), 
			       fitF[1]->GetParameter(2), fitF[1]->GetParError(2)));
   t3.DrawLatex(0.2,0.8,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF[2]->GetParameter(1), fitF[2]->GetParError(1), 
			       fitF[2]->GetParameter(2), fitF[2]->GetParError(2)));

   leg->Draw("same");
   chE->Print(Form((folder+"/energyPlots_layers%d_config%d.png").c_str(), nLayers, config), "png");
   chE->Print(Form((folder+"/energyPlots_layers%d_config%d.root").c_str(), nLayers, config), "root");
  */

   /////////////////////////////
   
   for(int iP=0; iP<5; ++iP)
     for(int iP=0; iP<5; ++iP){
       //       if(nLayers == 28 && iP >=4) continue;
       recoEnergy[iP]->Fit(Form("E%dR", energyP[iP]), "R");
     }
  
   TCanvas* chER = new TCanvas();
   chER->cd();
   recoEnergy[0]->GetXaxis()->SetTitle("#Sigma E_{i} / E_{beam}");
   recoEnergy[0]->GetXaxis()->SetRangeUser(10, 250.);
   recoEnergy[0]->GetYaxis()->SetRangeUser(0., 10.e+03);
   if(nLayers == 8){
     recoEnergy[0]->GetXaxis()->SetRangeUser(10., 250.);
     recoEnergy[0]->GetYaxis()->SetRangeUser(0., 3.e+03);
   }
   recoEnergy[0]->Draw();
   for(int iP=1; iP<5; ++iP){
     recoEnergy[iP]->Draw("same");
   }
   
   for(int iT=0; iT<5; ++iT){
     
     tg[0]->SetPoint(iT+1, energyP[iT], fitFRel[iT]->GetParameter(1)/energyP[iT]);
     tg[0]->SetPointError(iT+1, 0., fitFRel[iT]->GetParError(1) / energyP[iT]);

     tgL[0]->SetPoint(iT+1, energyP[iT], fitFRel[iT]->GetParameter(1));
     tgL[0]->SetPointError(iT+1, 0., sqrt(pow(fitFRel[iT]->GetParError(1), 2.)) );

     tgR[0]->SetPoint(iT+1, 1./energyP[iT], pow(fitFRel[iT]->GetParameter(2),2) / (1.* energyP[iT]) * 100.);
     tgR[0]->SetPointError(iT+1, 0., 100.*pow(fitFRel[iT]->GetParError(2),2.)/ (1.* energyP[iT]));

     tgS[0]->SetPoint(iT+1, energyP[iT], fitFRel[iT]->GetParameter(2) / fitFRel[iT]->GetParameter(1));
     tgS[0]->SetPointError(iT+1, 0., fitFRel[iT]->GetParameter(2) / fitFRel[iT]->GetParameter(1) * sqrt( pow(fitFRel[iT]->GetParError(2)/fitFRel[iT]->GetParameter(2), 2.) +
													 pow(fitFRel[iT]->GetParError(1)/fitFRel[iT]->GetParameter(1), 2.) ) );


     //std::cout << " >>> energyP[iT] = " << energyP[iT] << " fitF[iT]->GetParameter(1) / (1.* energyP[iT]) = " << fitF[iT]->GetParameter(1) / (1.* energyP[iT]) << std::endl;
     //     std::cout << " >>> energyP[iT] = " << energyP[iT] << " Y = " << pow(fitF[iT]->GetParameter(2),2) / (1.* energyP[iT]) * 100. << std::endl;

     tg[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgL[0]->SetPoint(tg[iT]->GetN(), 400, 300);
     tgR[0]->SetPoint(tg[iT]->GetN(), 400, 5);
     tgS[0]->SetPoint(tg[iT]->GetN(), 400, 5);
  }
   leg->Draw("same");
   if(nLayers != 28){
   chER->Print(Form((folder+"/energyPlots_layers%d_config%d.png").c_str(), nLayers, config), "png");
   chER->Print(Form((folder+"/energyPlots_layers%d_config%d.root").c_str(), nLayers, config), "root");
   }
   else{
   chER->Print(Form((folder+"/energyPlots_layers%d.png").c_str(), nLayers), "png");
   chER->Print(Form((folder+"/energyPlots_layers%d.root").c_str(), nLayers), "root");
   }



   /////////TGraph linearity                                                                                                                                                     
   TCanvas* tgLi = new TCanvas();
   tgLi->cd();
   tgL[0]->GetYaxis()->SetTitle("#Sigma E_{i} (GeV)");
   tgL[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
   tgL[0]->GetXaxis()->SetRangeUser(0, 300.);
   tg[0]->GetYaxis()->SetRangeUser(0., 1.2);
   tgL[0]->Draw("ap");
   if(nLayers != 28){
     tgLi->Print(Form((folder+"/energy_LinearityVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");
     tgLi->Print(Form((folder+"/energy_LinearityVsEnergy_layers%d_config%d.root").c_str(), nLayers, config), "root");
   }
   else{
     tgLi->Print(Form((folder+"/energy_LinearityVsEnergy_layers%d.png").c_str(), nLayers), "png");
     tgLi->Print(Form((folder+"/energy_LinearityVsEnergy_layers%d.root").c_str(), nLayers), "root");
   }


    /////////TGraph
    TCanvas* tgM = new TCanvas();
    tgM->cd();
    tg[0]->GetYaxis()->SetTitle("#Sigma E_{i} / E_{beam}");
    tg[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tg[0]->GetXaxis()->SetRangeUser(0, 300.);
    tg[0]->GetYaxis()->SetRangeUser(0., 1.01);
    tg[0]->GetYaxis()->SetRangeUser(0., 1.2);
    tg[0]->Draw("ap");
    if(nLayers != 28){
    tgM->Print(Form((folder+"/energy_MeanVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");
    tgM->Print(Form((folder+"/energy_MeanVsEnergy_layers%d_config%d.root").c_str(), nLayers, config), "root");
    }
    else{
      tgM->Print(Form((folder+"/energy_MeanVsEnergy_layers%d.png").c_str(), nLayers), "png");
      tgM->Print(Form((folder+"/energy_MeanVsEnergy_layers%d.root").c_str(), nLayers), "root");
    }

    /*
    /////////TGraph
    TCanvas* tgReso = new TCanvas();
    tgReso->cd();
    tgR[0]->GetYaxis()->SetTitle("#sigma(E) / E)^{2} (%)");
    tgR[0]->GetXaxis()->SetTitle("#gamma 1/E (GeV)");
    tgR[0]->GetXaxis()->SetRangeUser(0, 0.15);
    tgR[0]->GetYaxis()->SetRangeUser(0., 1.);
    tgR[0]->Draw("ap");
    tgR[1]->Draw("p, same");
    tgR[2]->Draw("p, same");
    legTGM->Draw("same");
    tgReso->Print((folder+"/energy_ResolutionVsEnergy.png").c_str(), "png");
    */


    TFile outSim(Form("SIM_resolution_layers%d_config%d.root", nLayers, config), "recreate");
    outSim.cd();
    tgS[0]->Write("resolution_GeV");
    tg[0]->Write("mean_GeV");
    tgL[0]->Write("linerity_GeV");
    outSim.Close();



    gStyle->SetOptFit(1);


    TF1* fitReso;
    //fit 3par all                                                                                                 
    //int type = 1;                                                                                                
    //fit 3par tail                                                                                                
    int type = 2;
    //fit 2par all                                                                                                 
    //int type = 3;                                                                                            
    if(type == 1){
      fitReso = new TF1("fitReso", "sqrt( pow([0]/sqrt(x), 2.) + pow([1]/x, 2.) + pow([2], 2.))", 10., 300.);
      fitReso->SetParName(0, "S");
      fitReso->SetParName(1, "N");
      fitReso->SetParName(2, "C");
      fitReso->SetParLimits(1, 0., 0.16);
      //    fitResoM->SetParLimits(2, 0.0001, 0.05);                                                               
    }
    if(type == 3){
      fitReso = new TF1("fitReso", "sqrt( pow([0]/sqrt(x), 2.) + pow([1], 2.))", 10., 300.);
      fitReso->SetParName(0, "S");
      fitReso->SetParName(1, "C");
      //    fitResoM->SetParLimits(1, 0., 0.05);                                                                  
      //    fitResoM->SetParLimits(1, 0.0001, 0.05);                                                              
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


    /////////TGraph   simple
    TCanvas* tgResoS = new TCanvas();
    tgResoS->cd();
    tgS[0]->GetYaxis()->SetTitle("#sigma( #Sigma E) / #Sigma E");
    tgS[0]->GetXaxis()->SetTitle("E_{beam} (GeV)");
    tgS[0]->GetXaxis()->SetRangeUser(0., 300.);
    tgS[0]->GetYaxis()->SetRangeUser(0., 0.3);
    tgS[0]->Draw("ap");
    fitReso->Draw("same");
    if(nLayers != 28){
      tgResoS->Print(Form((folder+"/energy_SResolutionVsEnergy_layers%d_config%d.png").c_str(), nLayers, config), "png");
      tgResoS->Print(Form((folder+"/energy_SResolutionVsEnergy_layers%d_config%d.root").c_str(), nLayers, config), "root");
    }
    else{
      tgResoS->Print(Form((folder+"/energy_SResolutionVsEnergy_layers%d.png").c_str(), nLayers), "png");
      tgResoS->Print(Form((folder+"/energy_SResolutionVsEnergy_layers%d.root").c_str(), nLayers), "root");
    }

    //    float X0value = 0;
    TGraphErrors* tM[5];
    for(int iT=0; iT<5; ++iT){
      tM[iT] = (TGraphErrors*)inF[iT]->Get("enLayer_MIP");
      tM[iT]->SetName(Form("enLayer_MIP_E%d", energyP[iT]));
      
      tM[iT]->SetMarkerStyle(20);
      tM[iT]->SetMarkerColor(iColors[iT]);
      tM[iT]->SetLineColor(iColors[iT]);
      tM[iT]->SetLineWidth(2);
      tM[iT]->SetLineStyle(iT);
    }

    TFile newOUT(Form("analyzed_SIM_layers%d_config%d.root", nLayers, config), "recreate");
    newOUT.cd();
    for(int iT=0; iT<5; ++iT){
      //      tM[iT]->SetName(Form("enLayer_MIP_E%d", energyP[iT]));
      tM[iT]->Write(Form("enLayer_MIP_E%d", energyP[iT]));
    }
    newOUT.Close();


    TCanvas* cTM = new TCanvas();
    cTM->cd();
    tM[0]->GetXaxis()->SetTitle("shower depth (X_{0})");
    tM[0]->GetYaxis()->SetTitle("#Sigma_{layer} E_{layer} (MIP)");
    tM[0]->GetXaxis()->SetRangeUser(0., 30.);
    tM[0]->GetYaxis()->SetRangeUser(0., 2.e3);
    tM[0]->Draw("ap");
    for(int iT=1; iT<5; ++iT){
      tM[iT]->Draw("p, same");
    }
    leg->Draw("same");
    if(nLayers != 28){
      cTM->Print(Form((folder+"/enLayer_MIP_nLayer%d_config%d.png").c_str(), nLayers, config), "png");
      cTM->Print(Form((folder+"/enLayer_MIP_nLayer%d_config%d.root").c_str(), nLayers, config), "root");
    }
    else{
      cTM->Print(Form((folder+"/enLayer_MIP_nLayer%d.png").c_str(), nLayers), "png");
      cTM->Print(Form((folder+"/enLayer_MIP_nLayer%d.root").c_str(), nLayers), "root");
        }
    return;
    ////////////////
    TFile* inF2[5];
    if(nLayers == 28 && config == 1){
      for(int iP=0; iP<4; ++iP)
	inF2[iP] = TFile::Open(Form("showers_ROOT/analyzed_28layers_Ele%d.root", energyP[iP]));
    }
    else if(nLayers == 8 && config == 2){
      for(int iP=0; iP<5; ++iP)
	inF2[iP] = TFile::Open(Form("showers_ROOT/analyzed_endShower_Ele%d.root", energyP[iP]));
    }
    else{
      for(int iP=0; iP<5; ++iP)
	inF2[iP] = TFile::Open(Form("showers_ROOT/analyzed_Ele%d.root", energyP[iP]));
    }


    std::cout << " >>> fatto " << std::endl;
    TGraphErrors* tEF[5];
    for(int iT=0; iT<5; ++iT){
      //      if(nLayers == 28 && iT >= 4) continue;
      tEF[iT] = (TGraphErrors*)inF2[iT]->Get("beamEfraction_layer");
      tEF[iT]->SetName(Form("beamEfractionvsL_E%d", energyP[iT]));
      
      tEF[iT]->SetMarkerStyle(20);
      tEF[iT]->SetMarkerColor(iColors[iT]);
      tEF[iT]->SetLineColor(iColors[iT]);
      tEF[iT]->SetLineWidth(2);
      tEF[iT]->SetLineStyle(iT);
    }
    


    TCanvas* tgBeamEF = new TCanvas();
    tgBeamEF->cd();
    tEF[0]->GetXaxis()->SetTitle("shower depth (X_{0})");
    tEF[0]->GetYaxis()->SetTitle("#Sigma_{layer} E_{layer} / E_{beam}");
    tEF[0]->GetXaxis()->SetRangeUser(0., 300.);
    tEF[0]->GetYaxis()->SetRangeUser(0., 1.1);
    tEF[0]->Draw("ap");
    for(int iT=1; iT<5; ++iT){
      //      if(nLayers == 28 && iT >= 4) continue;
      tEF[iT]->Draw("p, same");
    }
    legTGM->Draw("same");
    tgBeamEF->Print(Form((folder+"/beamEnergyFractionLayer%d_config%d.png").c_str(), nLayers, config), "png");
    

}
