#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include <TH1.h>
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

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit;

void simple_Data(int nLayer, int energyEle, int config = 1) {
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptTitle(0);

  float X0val[28];
  if(config == 1){
    float X0val_8[28] = {6.268, 1.131, 1.131, 1.362, 0.574, 1.301, 0.574, 2.42, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    for(int j=0; j<28; ++j){
      X0val[j] = X0val_8[j];
    }
  }
  else if(config == 2){
    float X0val_8[28] = {5.048, 3.412, 3.412, 2.866, 2.512, 1.625, 2.368, 6.021, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    for(int j=0; j<28; ++j){
      X0val[j] = X0val_8[j];
    }
  }


  TFile* inF;
  if(nLayer == 28) inF = TFile::Open(Form("energyLayer_ROOT/outF_28layers_Ele%d.root", energyEle));
  //  else if(config == 1) inF = TFile::Open(Form("~/public/perShilpi/cfg1_data_Mean/%dGeV.root", energyEle));
  //  else if(config == 2) inF = TFile::Open(Form("~/public/perShilpi/cfg2_data_Mean/%dGeV.root", energyEle));
  //  else if(config == 1) inF = TFile::Open(Form("~/public/perShilpi/cfg1_data_MPV/%dGeV.root", energyEle));
  //  else if(config == 2) inF = TFile::Open(Form("~/public/perShilpi/cfg2_data_MPV/%dGeV.root", energyEle));
  //bert
  //  else if(config == 1) inF = TFile::Open(Form("~/public/perShilpi/cfg1_data_MPV_BERT/%dGeV.root", energyEle));
  //else if(config == 2) inF = TFile::Open(Form("~/public/perShilpi/cfg2_data_MPV_BERT/%dGeV.root", energyEle));
  //+/- 2cm cut
  //  else if(config == 1) inF = TFile::Open(Form("~/public/perShilpi/cfg1_data_MPV_pm2cm/%dGeV.root", energyEle));
  //else if(config == 2) inF = TFile::Open(Form("~/public/perShilpi/cfg2_data_MPV_pm2cm/%dGeV.root", energyEle));
  else if(config == 1) inF = TFile::Open(Form("~/public/perShilpi/cfg1_data_MPV_logPos/%dGeV.root", energyEle));
  else if(config == 2) inF = TFile::Open(Form("~/public/perShilpi/cfg2_data_MPV_logPos/%dGeV.root", energyEle));


  TH1F* h_X_L1;
  TH1F* h_Y_L1;

  TH1F* enLayer[28];
  TH1F* enLayer_AbsW_Mip[28];
  TH1F* enLayer_AbsW_GeV[28];

  TH1F* enLayer_up[28];
  TH1F* enLayer_AbsW_Mip_up[28];
  TH1F* enLayer_AbsW_GeV_up[28];

  TH1F* enLayer_dw[28];
  TH1F* enLayer_AbsW_Mip_dw[28];
  TH1F* enLayer_AbsW_GeV_dw[28];

  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    if(i == 0){
      h_X_L1 = (TH1F*)inF->Get("LayerSumAnalyzer/X_L1");
      h_Y_L1 = (TH1F*)inF->Get("LayerSumAnalyzer/Y_L1");
    }

    enLayer[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d", i+1));
    enLayer_AbsW_Mip[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_AbsW_Mip", i+1));
    enLayer_AbsW_GeV[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_AbsW_GeV", i+1));

    enLayer_up[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_up", i+1));
    enLayer_AbsW_Mip_up[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_AbsW_Mip_up", i+1));
    enLayer_AbsW_GeV_up[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_AbsW_GeV_up", i+1));

    enLayer_dw[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_dw", i+1));
    enLayer_AbsW_Mip_dw[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_AbsW_Mip_dw", i+1));
    enLayer_AbsW_GeV_dw[i] = (TH1F*)inF->Get(Form("LayerSumAnalyzer/h_eAll_L%d_AbsW_GeV_dw", i+1));
  }

  std::cout << " >> presi histos " << std::endl;
  
  TGraphErrors* enLayer_MIP = new TGraphErrors();
  TGraphErrors* enLayer_MIPW = new TGraphErrors();
  TGraphErrors* enLayer_GeVW = new TGraphErrors();

  std::cout << " >> definiti graph " << std::endl;
  float X0value = 0;

  TGraphErrors* Mean_vsE = new TGraphErrors();
  TGraphErrors* Sigma_vsE = new TGraphErrors();
  //  TGraphErrors* beamEfraction_layer = new TGraphErrors();


  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;

    //enLayer[i]->Rebin(4);
    /*
    if(energyEle == 20) enLayer[i]->GetXaxis()->SetRangeUser(0., 0.02);
    if(energyEle == 100) enLayer[i]->GetXaxis()->SetRangeUser(0., 0.1);
    if(energyEle == 250) enLayer[i]->GetXaxis()->SetRangeUser(0., 0.2);
    if(energyEle == 200) enLayer[i]->GetXaxis()->SetRangeUser(0., 0.2);
    */

    float xMax = 0.;
    float xMaxVal = 0.;
    float xErr = 0.;

    float xMax_MW = 0.;
    float xMaxVal_MW = 0.;
    float xErr_MW = 0.;

    float xMax_GW = 0.;
    float xMaxVal_GW = 0.;
    float xErr_GW = 0.;

    for(int iB=1; iB<enLayer[i]->GetNbinsX()+1; ++iB){
      if(enLayer[i]->GetBinCenter(iB) == 0) continue;
      if(enLayer[i]->GetBinContent(iB) > xMaxVal){
	xMax = enLayer[i]->GetBinCenter(iB);
	xMaxVal = enLayer[i]->GetBinContent(iB);
	xErr = enLayer[i]->GetBinWidth(iB);
      }
      if(iB<enLayer_AbsW_Mip[i]->GetNbinsX()+1){
	if(enLayer_AbsW_Mip[i]->GetBinContent(iB) > xMaxVal_MW){
	  xMax_MW = enLayer_AbsW_Mip[i]->GetBinCenter(iB);
	  xMaxVal_MW = enLayer_AbsW_Mip[i]->GetBinContent(iB);
	  xErr_MW = enLayer_AbsW_Mip[i]->GetBinWidth(iB);
	}
      }
      if(iB<enLayer_AbsW_GeV[i]->GetNbinsX()+1){
	if(enLayer_AbsW_GeV[i]->GetBinContent(iB) > xMaxVal_GW){
	  xMax_GW = enLayer_AbsW_GeV[i]->GetBinCenter(iB);
	  xMaxVal_GW = enLayer_AbsW_GeV[i]->GetBinContent(iB);
	  xErr_GW = enLayer_AbsW_GeV[i]->GetBinWidth(iB);
	}
      }
    }

    //    std::cout << " >>> qui ok " << std::endl;

    //FIXME optimize fit
    /*

    //    if(i != 3) continue;
    // Declare observable x
    RooRealVar x("x", "x", 0., 0.02);
    if(energyEle == 20)  x.setMax(0.02);
    if(energyEle == 70) x.setMax(0.1);
    if(energyEle == 100) x.setMax(0.1);
    if(energyEle == 250) x.setMax(0.2);
    if(energyEle == 200) x.setMax(0.2);
    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist data("data", "data", x, Import(*enLayer[i]));
    */

    TF1* helpFit = new TF1("helpFit", "gaus", 0., 5.e+03);
    TF1* helpFit_MW = new TF1("helpFit_MW", "gaus", 0., 500.e+03);
    TF1* helpFit_GW = new TF1("helpFit_GW", "gaus", 0., 500.);
    enLayer[i]->Fit("helpFit", "RQ");
    enLayer_AbsW_Mip[i]->Fit("helpFit_MW", "RQ");
    enLayer_AbsW_GeV[i]->Fit("helpFit_GW", "RQ");


    /*
    // Construct observable
    RooRealVar ml("ml", "mean landau", helpFit->GetParameter(1), helpFit->GetParameter(1) * 0.8, enLayer[i]->GetMean());
    RooRealVar sl("sl", "sigma landau", helpFit->GetParameter(2)); //, helpFit->GetParameter(2)*0.5, helpFit->GetParameter(2));
    RooLandau landau("lx", "lx", x, ml, sl);
    //  landau.fitTo(data);  
    //    std::cout << " L = " << ml.getVal() << " +/- " << ml.getError() << std::endl;

    // Construct gauss(t,mg,sg)
    RooRealVar mg("mg", "mg", 0.);
    //RooRealVar mg("mg", "mg",  enLayer[i]->GetMean(), 0., 0.2);
    RooRealVar sg("sg", "sg", enLayer[i]->GetRMS(), enLayer[i]->GetRMS()/2., enLayer[i]->GetRMS()*2.);
    //    RooRealVar sg("sg", "sg", enLayer[i]->GetRMS());
    RooGaussian gauss("gauss", "gauss", x, mg, sg);
    
  
    // C o n s t r u c t   c o n v o l u t i o n   p d f 
    // ---------------------------------------
    // Set #bins to be used for FFT sampling to 10000
    x.setBins(100000,"cache"); 

    // Construct landau (x) gauss
    RooFFTConvPdf lxg("lxg", "landau (X) gauss", x, landau, gauss);
    // Fit gxlx to data
    //gauss.fitTo(*data);
    RooFitResult* info = lxg.fitTo(data);
    //RooFitResult* info = landau.fitTo(data);
    // Plot data, landau pdf, landau (X) gauss pdf
    RooPlot* frame1 = x.frame(Title("landau (x) gauss convolution"));
    data.plotOn(frame1);
    lxg.plotOn(frame1);
    //landau.plotOn(frame1, LineStyle(kDashed));
    //gauss.plotOn(frame1, LineColor(kRed));
         

    TF1 *f = lxg.asTF( RooArgList(x) );
    float convMax = f->GetMaximumX(); 
     
    
    //    info->Print();
    std::cout << " convMax = " << convMax << std::endl;
    std::cout << " max = " << xMax << std::endl;
    //    std::cout << " Gsig = " << sg.getVal() << " +/- " << sg.getError() << std::endl;
    */

    // Draw frame on canvas
    TCanvas* c2 = new TCanvas("fit", "", 600, 600);
    c2->cd();
    enLayer[i]->GetXaxis()->SetRangeUser(0., 3.e+03);
    enLayer[i]->Draw();
    c2->Print(Form("fitFolder/histoLayer%d_energy%d.png", i, energyEle), ".png");


    X0value += X0val[i];
    /*
    enLayer_MIP->SetPoint(i+1, X0value, helpFit->GetParameter(1));
    //    enLayer_MIP->SetPointError(i+1, 0, helpFit->GetParError(1));
    enLayer_MIPW->SetPoint(i+1, X0value, helpFit_MW->GetParameter(1));
    enLayer_GeVW->SetPoint(i+1, X0value, helpFit_GW->GetParameter(1));
    */

    std::cout << " >> riempio graph " << i << " val " << enLayer[i]->GetMean() 
	      << " val up " << enLayer_up[i]->GetMean() 
	      << " val dw " << enLayer_dw[i]->GetMean() 
	      << std::endl;

    enLayer_MIP->SetPoint(i+1, X0value, enLayer[i]->GetMean());      
    enLayer_MIP->SetPointError(i+1, 0, std::abs(enLayer_up[i]->GetMean() - enLayer_dw[i]->GetMean())/2.);      
    //    enLayer_MIP->SetPointError(i+1, 0, helpFit->GetParError(1));  


    std::cout << " >> riempio graph " << i << " val " << enLayer_AbsW_Mip[i]->GetMean() 
	      << " val up " << enLayer_AbsW_Mip_up[i]->GetMean() 
	      << " val dw " << enLayer_AbsW_Mip_dw[i]->GetMean() 
	      << std::endl;

    enLayer_MIPW->SetPoint(i+1, X0value, enLayer_AbsW_Mip[i]->GetMean());  
    enLayer_MIPW->SetPointError(i+1, 0, std::abs(enLayer_AbsW_Mip_up[i]->GetMean() - enLayer_AbsW_Mip_dw[i]->GetMean())/2.);  

    enLayer_GeVW->SetPoint(i+1, X0value, enLayer_AbsW_GeV[i]->GetMean());  
    enLayer_GeVW->SetPointError(i+1, 0, std::abs(enLayer_AbsW_GeV_up[i]->GetMean() - enLayer_AbsW_GeV_dw[i]->GetMean()) /2.);  
  }

  std::cout << " fatto " << std::endl;

  TH1F* sumAll_Layer = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all");
  TH1F* sumAll_Layer_MIP = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_AbsW_Mip");
  TH1F* sumAll_Layer_GeV = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_AbsW_GeV");

  TH1F* sumAll_Layer_up = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_up");
  TH1F* sumAll_Layer_MIP_up = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_AbsW_Mip_up");
  TH1F* sumAll_Layer_GeV_up = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_AbsW_GeV_up");

  TH1F* sumAll_Layer_dw = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_dw");
  TH1F* sumAll_Layer_MIP_dw = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_AbsW_Mip_dw");
  TH1F* sumAll_Layer_GeV_dw = (TH1F*)inF->Get("LayerSumAnalyzer/h_eAll_all_AbsW_GeV_dw");


  // TGraph* Mean_vsE = new TGraph();
  // TGraph* Sigma_vsE = new TGraph();


  h_X_L1->Rebin(4);
  h_Y_L1->Rebin(4);

  TF1* fitF = new TF1("fitF", "gaus", -6, 6);
  h_X_L1->Fit("fitF");

  TCanvas* cX = new TCanvas();
  cX->cd();
  h_X_L1->GetXaxis()->SetTitle("X");
  h_X_L1->Draw();
  fitF->Draw("same");
  cX->Print(Form("position/h_X_L1_Ele%d_cfg%d.png", energyEle, config), "png");

  h_Y_L1->Fit("fitF");
  TCanvas* cY = new TCanvas();
  cY->cd();
  h_Y_L1->GetXaxis()->SetTitle("Y");
  h_Y_L1->Draw();
  fitF->Draw("same");
  cY->Print(Form("position/h_Y_L1_Ele%d_cfg%d.png", energyEle, config), "png");



  TFile* outF;
  if(nLayer == 28) outF = new TFile(Form("showers_ROOT_data/analyzed_28layers_Ele%d.root", energyEle), "recreate");
  else if(config == 1) outF = new TFile(Form("showers_ROOT_data/analyzed_Ele%d.root", energyEle), "recreate");
  else if(config == 2) outF = new TFile(Form("showers_ROOT_data/analyzed_endShower_Ele%d.root", energyEle), "recreate");

  outF->cd();
  sumAll_Layer->Write();
  sumAll_Layer_MIP->Write();
  sumAll_Layer_GeV->Write();

  sumAll_Layer_up->Write();
  sumAll_Layer_MIP_up->Write();
  sumAll_Layer_GeV_up->Write();

  sumAll_Layer_dw->Write();
  sumAll_Layer_MIP_dw->Write();
  sumAll_Layer_GeV_dw->Write();

  enLayer_MIP->Write("enLayer_MIP");
  enLayer_MIPW->Write("enLayer_MIPW");
  enLayer_GeVW->Write("enLayer_GeVW");
  outF->Close();

}
  
