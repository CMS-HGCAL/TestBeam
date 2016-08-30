#include <memory>
#include <iostream>
#include "TStyle.h"
#include "TH2Poly.h"
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


void DumpPlotsReco(TString inputFileName, TString outputFolder, Int_t runNumber) {
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->Range(0, 0, 1, 1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);

  
  // char filename[50], dirname[100], histname[50];
  // int run = 2;
  // int spill = 3;
  const int nLayers = 1;
  const int sensorsize = 128;
  const int Sensor_Iu = 0;
  const int Sensor_Iv = 0;

  HGCalTBTopology IsCellValid;

  // for(int iii=0; iii<=run; iii++){
  //   for(int jjj=1; jjj<=spill;jjj++){   
  // if(iii != 0) sprintf(filename,"test_DigiAndRechitPlotter_TB_Run%i_Spill%i.root",iii,jjj);
  // else sprintf(filename,"test_DigiAndRechitPlotter_TB_Cumulative.root");
  TFile *inputFile = TFile::Open(inputFileName);
  inputFile->cd("hgcaltbrechitsplotter");

  TString outputPlotName, histogramName;

  for(int layerCounter=1; layerCounter<= nLayers; layerCounter++){
    histogramName = Form("FullLayer_RecHits_Layer%i",layerCounter);
    // sprintf(histname,"FullLayer_RecHits_Layer%i",layerCounter);
    TH2Poly* Overview_Rechits_Layer  = (TH2Poly*) inputFile->FindObjectAny(histogramName);
    Overview_Rechits_Layer->GetYaxis()->SetTitle("Y[cm]");
    Overview_Rechits_Layer->GetYaxis()->CenterTitle(1);
    Overview_Rechits_Layer->Draw("colztext");
    ReverseXAxis(Overview_Rechits_Layer);
    // if(iii != 0) sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_RecHits_Layer%i.png",iii,jjj,layerCounter);
    // else sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_RecHits_Layer%i.png",layerCounter);
    outputPlotName = (outputFolder + Form("/Overview/")) + histogramName + Form("_%06d.png", runNumber);
    c1->SaveAs(outputPlotName);
    delete Overview_Rechits_Layer;

    histogramName = Form("FullLayer_Occupancy_Layer%i", layerCounter);
    TH2Poly* Overview_Occupancy_Layer  = (TH2Poly*) inputFile->FindObjectAny(histogramName);
    Overview_Occupancy_Layer->GetYaxis()->SetTitle("Y[cm]");
    Overview_Occupancy_Layer->GetYaxis()->CenterTitle(1);
    Overview_Occupancy_Layer->Draw("colztext");
    ReverseXAxis(Overview_Occupancy_Layer);
    // if(iii != 0) sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_Occupancy_Layer%i.png",iii,jjj,layerCounter);
    // else sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_Occupancy_Layer%i.png",layerCounter);
    outputPlotName = (outputFolder + Form("/Overview/")) + histogramName + Form("_%06d.png", runNumber);
    c1->SaveAs(outputPlotName);
    delete Overview_Occupancy_Layer;
    
    histogramName = Form("FullLayer_RecHits_Layer%i_Summed",layerCounter);
    TH1F* Overview_Rechits_Layer_Summed  = (TH1F*) inputFile->FindObjectAny(histogramName);
    Overview_Rechits_Layer_Summed->Draw();
    // if(iii != 0) sprintf(outputPlotName,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_RecHits_Layer%i_Summed.png",iii,jjj,layerCounter);
    // else sprintf(outputPlotName,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_RecHits_Layer%i_Summed.png",layerCounter);
    outputPlotName = (outputFolder + Form("/Overview/")) + histogramName + Form("_%06d.png", runNumber);
    c1->SaveAs(outputPlotName);
    delete Overview_Rechits_Layer_Summed;
	
    for(int iv = -7; iv < 8; iv++){
      for(int iu = -7; iu < 8; iu++){
        if(!IsCellValid.iu_iv_valid(nLayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
        histogramName = Form("Cell_RecHits_u_%i_v_%i_Layer%i",iu,iv,layerCounter);
        TH1F* Cell_Hist_U_V  = (TH1F*) inputFile->FindObjectAny(histogramName);
        // if(iii != 0) sprintf(outputPlotName,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Detailed/Cell_RecHits_u_%i_v_%i_Layer_%i.png",iii,jjj,iu,iv,layerCounter);
        // else sprintf(outputPlotName,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Detailed/Cell_RecHits_u_%i_v_%i_Layer_%i.png",iu,iv,layerCounter);
        outputPlotName = (outputFolder + Form("/Detailed/")) + histogramName + Form("_%06d.png", runNumber);
        Cell_Hist_U_V->Draw();
        c1->SaveAs(outputPlotName);
        delete Cell_Hist_U_V;
      }//loop over iu ends here
    }//loop over iv ends here

  }//loop over layer ends here
  //   }//loop over spill ends here
  // }//loop over run ends here
}
