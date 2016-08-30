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

void DumpPlotsDigi(TString inputFileName, TString outputFolder, Int_t runNumber) {
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->Range(0, 0, 1, 1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  TFile *inputFile = TFile::Open(inputFileName);
  inputFile->cd("hgcaltbdigisplotter");
  TString outputPlotFullPath, plotName;
  const int nLayers = 1;
  const int nSkirocs = 2;
  const int nChannelsPerSkiroc = 64;
      
  for (int gainCounter = 0; gainCounter <= 1; ++gainCounter) {
    for (int layerCounter = 1; layerCounter <= nLayers; ++layerCounter) {
      plotName = Form("FullLayer_ADC%i_Layer%i", gainCounter, layerCounter);
      TH2Poly* sumADCCounts = (TH2Poly*) inputFile->FindObjectAny(plotName);
      sumADCCounts->GetYaxis()->SetTitle("Y[cm]");
      sumADCCounts->GetYaxis()->CenterTitle(1);
      sumADCCounts->Draw("colztext");
      ReverseXAxis(sumADCCounts);
      // if(iii != 0) sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_RecHits_Layer%i.png",iii,jjj,layerCounter);
      // else sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_RecHits_Layer%i.png",layerCounter);
      outputPlotFullPath = (outputFolder + Form("/Overview/")) + plotName + Form("_%06d.png", runNumber);
      c1->SaveAs(outputPlotFullPath);
      delete sumADCCounts;

      plotName = Form("FullLayer_ADC%i_Layer%i_profile", gainCounter, layerCounter);
      TProfile* allADCCounts = (TProfile*) inputFile->FindObjectAny(plotName);
      allADCCounts->GetXaxis()->SetTitle("Channel number");
      allADCCounts->GetYaxis()->SetTitle("ADC Counts");
      allADCCounts->Draw();
      outputPlotFullPath = (outputFolder + Form("/Overview/")) + plotName + Form("_%06d.png", runNumber);
      c1->SaveAs(outputPlotFullPath);
      delete allADCCounts;
    } // ends loop over layers
    for (int skirocCounter = 1; skirocCounter <= nSkirocs; ++skirocCounter) {
      for (int channelCounter = 0; channelCounter < nChannelsPerSkiroc; ++channelCounter) {
        plotName = Form("Ski_%i_Channel_%i_ADC%i", skirocCounter, channelCounter, gainCounter);
        TH1F* digisDistribution = (TH1F*) inputFile->FindObjectAny(plotName);
        digisDistribution->GetXaxis()->SetTitle("ADC Counts");
        digisDistribution->Draw();
        outputPlotFullPath = (outputFolder + Form("/Detailed/")) + plotName + Form("_%06d.png", runNumber);
        c1->SaveAs(outputPlotFullPath);
        delete digisDistribution;
      } // ends loop over channels in given skiroc
    } // ends loop over all skirocs
  } // ends loop over low/high gain
}
