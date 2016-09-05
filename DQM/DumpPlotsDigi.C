#include <memory>
#include <iostream>
#include "TStyle.h"
#include "TProfile.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGaxis.h"
//#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
//#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/src/HGCalTBCellVertices.cc"
#include "HGCal/Geometry/src/HGCalTBTopology.cc"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

void ReverseXAxis (TH1 *h)
{
  h->GetXaxis()->SetLabelOffset(999);
  h->GetXaxis()->SetTickLength(0);

  gPad->Update();
  TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
                               gPad->GetUymin(),
                               gPad->GetUxmin(),
                               gPad->GetUymin(),
                               h->GetXaxis()->GetXmin(),
                               h->GetXaxis()->GetXmax(),
                               510,"-");
  newaxis->SetLabelOffset(-0.03);
  newaxis->SetTitle("X[cm]");
  newaxis->CenterTitle(1);
  newaxis->Draw();
}

void histo_style(TH1F *h1, double xlow, double xhi, int rebin =1, bool changerange=true)
{
  h1->Rebin(rebin);
  if(changerange)
    h1->GetXaxis()->SetRangeUser(xlow, xhi);
}

void DumpPlotsDigi(TString inputFileName, TString outputFolder, Int_t runNumber, Int_t nSpills) {
  
  std::cout << "DumpPlotsDigi called for inputFileName: " << inputFileName << ", outputFolder: " << outputFolder << ", runNumber: " << runNumber << ", nSpills: " << nSpills << std::endl;

  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);

  TCanvas *plotCanvas = new TCanvas("plotCanvas", "plotCanvas", 600, 600);
  plotCanvas->Range(0, 0, 1, 1);
  plotCanvas->SetFillColor(0);
  plotCanvas->SetBorderMode(0);
  plotCanvas->SetBorderSize(2);
  plotCanvas->SetFrameBorderMode(0);

  const int nLayers = 8;
  int nSkirocs = nLayers*2;
  const int nChannelsPerSkiroc = 64;

  int nEvents = nSpills*EVENTSPERSPILL*nLayers;

  TFile *inputFile = TFile::Open(inputFileName);
  inputFile->cd("hgcaltbdigisplotter");
  TString outputPlotFullPath, plotName;
        
  for (int gainCounter = 0; gainCounter <= 1; ++gainCounter) {
    for (int layerCounter = 1; layerCounter <= nLayers; ++layerCounter) {
      plotName = Form("FullLayer_ADC%i_Layer%i", gainCounter, layerCounter);
      TH2Poly* sumADCCounts = (TH2Poly*) inputFile->FindObjectAny(plotName);
      if (not(NULL == sumADCCounts)) {
        sumADCCounts->GetYaxis()->SetTitle("Y[cm]");
        sumADCCounts->GetYaxis()->CenterTitle(1);
        sumADCCounts->Draw("colztext");
        ReverseXAxis(sumADCCounts);
        // if(iii != 0) sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_RecHits_Layer%i.png",iii,jjj,layerCounter);
        // else sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_RecHits_Layer%i.png",layerCounter);
        outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
        plotCanvas->SaveAs(outputPlotFullPath);
      }
      delete sumADCCounts;

      plotName = Form("FullLayer_ADC%i_Layer%i_profile", gainCounter, layerCounter);
      TProfile* allADCCounts = (TProfile*) inputFile->FindObjectAny(plotName);
      if (not(NULL == allADCCounts)) {
        allADCCounts->GetXaxis()->SetTitle("Channel number");
        allADCCounts->GetYaxis()->SetTitle("ADC Counts");
        allADCCounts->Draw();
        outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
        plotCanvas->SaveAs(outputPlotFullPath);
        delete allADCCounts;
      }
    } // ends loop over layers
    for (int skirocCounter = 1; skirocCounter <= nSkirocs; ++skirocCounter) {
      for (int channelCounter = 0; channelCounter < nChannelsPerSkiroc; ++channelCounter) {
        plotName = Form("Ski_%i_Channel_%i_ADC%i", skirocCounter, channelCounter, gainCounter);
        TH1F* digisDistribution = (TH1F*) inputFile->FindObjectAny(plotName);
        digisDistribution->GetXaxis()->SetTitle("ADC Counts");
        digisDistribution->Draw();
        outputPlotFullPath = (outputFolder + Form("/Detailed/")) + plotName + Form("_%06d.png", runNumber);
        plotCanvas->SaveAs(outputPlotFullPath);
        delete digisDistribution;
      } // ends loop over channels in given skiroc
    } // ends loop over all skirocs
  } // ends loop over low/high gain

  delete inputFile;

  delete plotCanvas;
}
