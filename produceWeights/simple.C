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

void simple(int nLayer, int energyEle, int config) {
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
      float dEdX_weights_8[28] = {33.151, 13.261, 14.247, 9.885, 9.844, 9.844, 16.416, 28.296, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
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


  // float weights2GeV = 1.e-03;
  // //weights2GeV = weights2GeV / 0.828;
  // //  float MIP2GeV_sim = 77.52e-06;
  // float MIP2GeV_sim = 62.67e-06;


  float MIP2GeV_sim = 5.28e-05;
  float weights2GeV = 1.e-03;
  float weights2MIP = 52.8 / 63.6;
      //  float MIP2GeV_sim = 64.19e-06;         
      //  float MIP2GeV_sim = 53.16e-06;         
      //  float MIP2GeV_sim = 53.16e-06;         
      //  weights2GeV = weights2GeV / 0.828;     
  float G4Escale = 1.;  //0.83;   // >>> fix scale i

  TFile* inF;
  if(nLayer == 28) inF = TFile::Open(Form("energyLayer_ROOT/outF_28layers_Ele%d.root", energyEle));
  else if(config == 1) inF = TFile::Open(Form("energyLayer_ROOT/outF_Ele%d.root", energyEle));
  else if(config == 2) inF = TFile::Open(Form("energyLayer_ROOT/outF_endShower_Ele%d.root", energyEle));

  TH1F* enLayer[28];
  for(int i=0; i<28; ++i){
    if(nLayer != 28 && i >= 8) continue;
    enLayer[i] = (TH1F*)inF->Get(Form("enLayer%d", i+1));
  }


  TGraph* enLayer_MPV = new TGraph();
  TGraph* enLayer_Mean = new TGraph();
  TGraph* enLayer_MIP = new TGraph();
  TGraph* enLayer_MeanW = new TGraph();
  TGraph* enLayer_MeanWMIP = new TGraph();
  TGraphErrors* beamEfraction_layer = new TGraphErrors();

  float totEnergyErr = 0;
  float totEnergy = 0;
  float totMIP = 0;

  float X0value = 0;

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

    for(int iB=1; iB<enLayer[i]->GetNbinsX()+1; ++iB){
      if(enLayer[i]->GetBinContent(iB) > xMaxVal){
	xMax = enLayer[i]->GetBinCenter(iB);
	xMaxVal = enLayer[i]->GetBinContent(iB);
	xErr = enLayer[i]->GetBinWidth(iB);
      }
    }

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


    TF1* helpFit = new TF1("helpFit", "landau", 0., 0.2);
    enLayer[i]->Fit("helpFit");


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
     
    // Draw frame on canvas
    TCanvas* c2 = new TCanvas("fit", "", 600, 600);
    gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
    c2->Print(Form("fitFolder/histoLayer%d_energy%d.png", i, energyEle), ".png");

    //    info->Print();
    std::cout << " convMax = " << convMax << std::endl;
    std::cout << " max = " << xMax << std::endl;
    //    std::cout << " Gsig = " << sg.getVal() << " +/- " << sg.getError() << std::endl;
    
    //    float estimateE = ml.getVal();
    float estimateE = xMax;
    //    float estimateEerr = min(abs(f1->GetParameter(1) - f2->GetParameter(1)), abs(f1->GetParameter(1) - f3->GetParameter(1)) );
    //    float estimateEerr = ml.getError();
    //float estimateEerr = sqrt(pow(ml.getError(), 2) + pow(estimateE - xMax, 2));
    //    float estimateEerr = sqrt(pow(ml.getError(), 2) + pow(estimateE - ml.getVal(), 2));
    //    float estimateEerr = min(enLayer[i]->GetRMS(), abs(estimateE - ml.getVal()));
    //    float estimateEerr = enLayer[i]->GetRMS();
    float estimateEerr = TMath::Max(xErr, xMax - convMax);

    std::cout << " xMax = " << xMax << " +/- " << estimateEerr << std::endl;

    Double_t chi2 = frame1->chiSquare("lxg", "data", 3);
    std::cout << " chi2 lxg = " << chi2 << std::endl;


    X0value += X0val[i];

    enLayer_MPV->SetPoint(i+1, X0value, estimateE);
    enLayer_Mean->SetPoint(i+1, X0value, estimateE);
    enLayer_MIP->SetPoint(i+1, X0value, estimateE/MIP2GeV_sim);
    
    enLayer_MeanW->SetPoint(i+1, X0value, estimateE / G4Escale * ( 1./MIP2GeV_sim * weights2GeV * weights2MIP * dEdX_weights[i] + 1.));
    //enLayer_MeanW->SetPoint(i+1, X0value, estimateE ); // * ( 1./MIP2GeV_sim * weights2GeV * dEdX_weights[i] + 1.));
    //    enLayer_MeanWMIP->SetPoint(i+1, X0value, estimateE * ( 1./MIP2GeV_sim * weights2GeV * dEdX_weights[i] + 1.) / MIP2GeV_sim );
    //    enLayer_MeanWMIP->SetPoint(i+1, X0value, estimateE / MIP2GeV_sim );

    /*
    totEnergy += estimateE * ( 1./MIP2GeV_sim * weights2GeV * dEdX_weights[i] + 1.);
    totMIP += estimateE * ( 1./MIP2GeV_sim * weights2GeV * dEdX_weights[i] + 1.) / MIP2GeV_sim;
    */
    //    totEnergy += estimateE;
    
    totMIP += estimateE / MIP2GeV_sim;

    //    totEnergyErr =  sqrt(pow(totEnergyErr, 2) + pow(estimateEerr * ( 1./MIP2GeV_sim * weights2GeV * dEdX_weights[i] + 1.), 2));
    totEnergy += estimateE / G4Escale * ( 1./MIP2GeV_sim * weights2GeV * weights2MIP * dEdX_weights[i] + 1.);
    totEnergyErr =  sqrt( pow(estimateEerr / G4Escale * ( 1./MIP2GeV_sim * weights2GeV * weights2MIP * dEdX_weights[i] + 1.), 2));
    // totEnergyErr =  sqrt( pow(estimateEerr, 2));
    beamEfraction_layer->SetPoint(i+1, X0value, totEnergy / (1.*energyEle));
    beamEfraction_layer->SetPointError(i+1, 0., totEnergyErr);
  }


  TFile* outF;
  if(nLayer == 28) outF = new TFile(Form("showers_ROOT/analyzed_28layers_Ele%d.root", energyEle), "recreate");
  else if(config == 1) outF = new TFile(Form("showers_ROOT/analyzed_Ele%d.root", energyEle), "recreate");
  else if(config == 2) outF = new TFile(Form("showers_ROOT/analyzed_endShower_Ele%d.root", energyEle), "recreate");

  outF->cd();
  enLayer_MPV->Write("enLayer_MPV");
  enLayer_Mean->Write("enLayer_Mean");
  enLayer_MIP->Write("enLayer_MIP");
  enLayer_MeanW->Write("enLayer_MeanW");
  enLayer_MeanWMIP->Write("enLayer_MeanWMIP");
  beamEfraction_layer->Write("beamEfraction_layer");
  outF->Close();


  std::cout << " >>> beam energy in MIP = " << 1.* energyEle / MIP2GeV_sim << " in energy = " << 1.* energyEle << std::endl;
  std::cout << " >>> " << nLayer << " layers >> sum MIP = " << totMIP << " sum energy = " << totEnergy << std::endl;
  std::cout << " >>> " << nLayer << " layers >> sum/beam = " << totEnergy/(1.* energyEle) << " energy and in MIP = " << totMIP / (1.* energyEle / MIP2GeV_sim) << std::endl;

}
  
