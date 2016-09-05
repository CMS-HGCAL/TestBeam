#include <memory>
#include <iostream>
#include "TStyle.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TH2F.h"
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

void DumpPlotsReco(TString inputFileName, TString outputFolder, Int_t runNumber, Int_t nSpills) {

  std::cout << "DumpPlotsReco called for inputFileName: " << inputFileName << ", outputFolder: " << outputFolder << ", runNumber: " << runNumber << ", nSpills: " << nSpills << std::endl;

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
  // const int nChannels = 64;

  int nEvents = nSpills*EVENTSPERSPILL*nLayers;

  TFile *inputFile = TFile::Open(inputFileName);
  inputFile->cd("hgcaltbrechitsplotter_highgain_correlation_cm");
  TString outputPlotFullPath, plotName;

  plotName = Form("Noise_2D_Profile_Layer");
  TH2F* h_Noise_2D_Profile_Layer  = (TH2F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_Noise_2D_Profile_Layer)) {
    h_Noise_2D_Profile_Layer->Draw("colz");
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_Noise_2D_Profile_Layer;

  plotName = Form("AllCells_Ped");
  TH1F* h_AllCells_Ped  = (TH1F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_AllCells_Ped)) {
    h_AllCells_Ped->Draw();
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_AllCells_Ped;

  plotName = Form("AllCells_CM");
  TH1F* h_AllCells_CM  = (TH1F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_AllCells_CM)) {
    h_AllCells_CM->Draw();
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_AllCells_CM;

  for(int layerCounter=0; layerCounter < nLayers; ++layerCounter){
    plotName = Form("Full_Cell_Layer_%d", layerCounter);
    TH1F* h_Full_Cell  = (TH1F*) inputFile->FindObjectAny(plotName);
    if (not(NULL == h_Full_Cell)) {
      h_Full_Cell->Draw();
      outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
      plotCanvas->SaveAs(outputPlotFullPath);
    }
    delete h_Full_Cell;

    plotName = Form("Half_Cell_Layer_%d", layerCounter);
    TH1F* h_Half_Cell  = (TH1F*) inputFile->FindObjectAny(plotName);
    if (not(NULL == h_Half_Cell)) {
      h_Half_Cell->Draw();
      outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
      plotCanvas->SaveAs(outputPlotFullPath);
    }
    delete h_Half_Cell;

    plotName = Form("MB_Cell_Layer_%d", layerCounter);
    TH1F* h_MB_Cell  = (TH1F*) inputFile->FindObjectAny(plotName);
    if (not(NULL == h_MB_Cell)) {
      h_MB_Cell->Draw();
      outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
      plotCanvas->SaveAs(outputPlotFullPath);
    }
    delete h_MB_Cell;

    plotName = Form("Calib_Pads_Layer_%d", layerCounter);
    TH1F* h_Calib_Pads  = (TH1F*) inputFile->FindObjectAny(plotName);
    if (not(NULL == h_Calib_Pads)) {
      h_Calib_Pads->Draw();
      outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
      plotCanvas->SaveAs(outputPlotFullPath);
    }
    delete h_Calib_Pads;

    plotName = Form("Merged_Cell_Layer_%d", layerCounter);
    TH1F* h_Merged_Cell  = (TH1F*) inputFile->FindObjectAny(plotName);
    if (not(NULL == h_Merged_Cell)) {
      h_Merged_Cell->Draw();
      outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
      plotCanvas->SaveAs(outputPlotFullPath);
    }
    delete h_Merged_Cell;

    // for (int skirocCounter = 1; skirocCounter <= nSkirocs; ++skirocCounter) {
    //   for (int channelCounter = 0; channelCounter <= nChannels; ++channelCounter) {
    //     plotName = Form("Ski_%d_Channel_%d_Layer_%d", skirocCounter, channelCounter, layerCounter);
    //     TH1F* h_Ski  = (TH1F*) inputFile->FindObjectAny(plotName);
    //     if (not(NULL == h_Ski)) {
    //       h_Ski->Draw();
    //       outputPlotFullPath = (outputFolder + Form("/Detailed/")) + plotName + Form("_%06d.png", runNumber);
    //       plotCanvas->SaveAs(outputPlotFullPath);
    //     }
    //     delete h_Ski;
    //   } // ends loop over channels
    // } // ends loop over skirocs

  } //loop over layers ends here

  inputFile->cd("hgcaltbrechitsplotter_highgain_new");

  plotName = Form("Covar_hist");
  TH2F* h_Covar_hist  = (TH2F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_Covar_hist)) {
    h_Covar_hist->Draw("colz");
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_Covar_hist;

  plotName = Form("Correl_hist");
  TH2F* h_Correl_hist  = (TH2F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_Correl_hist)) {
    h_Correl_hist->Draw("colz");
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_Correl_hist;

  plotName = Form("DiffIJ_hist");
  TH2F* h_DiffIJ_hist  = (TH2F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_DiffIJ_hist)) {
    h_DiffIJ_hist->Draw("colz");
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_DiffIJ_hist;

  plotName = Form("CG_X");
  TH1F* h_CG_X  = (TH1F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_CG_X)) {
    h_CG_X->Draw();
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_CG_X;

  plotName = Form("CG_Y");
  TH1F* h_CG_Y  = (TH1F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_CG_Y)) {
    h_CG_Y->Draw();
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_CG_Y;

  // plotName = Form("AllCells_Ped");
  // TH1F* h_AllCells_Ped  = (TH1F*) inputFile->FindObjectAny(plotName);
  // if (not(NULL == h_AllCells_Ped)) {
  //   h_AllCells_Ped->Draw();
  //   outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
  //   plotCanvas->SaveAs(outputPlotFullPath);
  // }
  // delete h_AllCells_Ped;

  // plotName = Form("AllCells_CM");
  // TH1F* h_AllCells_CM  = (TH1F*) inputFile->FindObjectAny(plotName);
  // if (not(NULL == h_AllCells_CM)) {
  //   h_AllCells_CM->Draw();
  //   outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
  //   plotCanvas->SaveAs(outputPlotFullPath);
  // }
  // delete h_AllCells_CM;

  plotName = Form("Sum_Cluster_ADC");
  TH1F* h_Sum_Cluster_ADC  = (TH1F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_Sum_Cluster_ADC)) {
    h_Sum_Cluster_ADC->Draw();
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_Sum_Cluster_ADC;

  plotName = Form("Sum_Cluster_Max");
  TH1F* h_Sum_Cluster_Max  = (TH1F*) inputFile->FindObjectAny(plotName);
  if (not(NULL == h_Sum_Cluster_Max)) {
    h_Sum_Cluster_Max->Draw();
    outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    plotCanvas->SaveAs(outputPlotFullPath);
  }
  delete h_Sum_Cluster_Max;

  for(int layerCounter=1; layerCounter <= nLayers; ++layerCounter){
    plotName = Form("FullLayer_RecHits_Layer%d_Summed", layerCounter);
    TH1F* h_FullLayer_RecHits  = (TH1F*) inputFile->FindObjectAny(plotName);
    if (not(NULL == h_FullLayer_RecHits)) {
      h_FullLayer_RecHits->Draw();
      outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
      plotCanvas->SaveAs(outputPlotFullPath);
    }
    delete h_FullLayer_RecHits;
  } //loop over layers ends here

  for (int eventCenturyCounter = 0; eventCenturyCounter < int(nEvents/100); ++eventCenturyCounter) {
    int eventNumber = eventCenturyCounter*100;
    int layerNumber = 1+((int(eventNumber/EVENTSPERSPILL))%nLayers);
    plotName = Form("FullLayer_ADC0_Layer%d_Event%d", layerNumber, eventNumber);
    TH2Poly* h_EventDisplay  = (TH2Poly*) inputFile->FindObjectAny(plotName);
    if (not(NULL == h_EventDisplay)) {
      h_EventDisplay->GetYaxis()->SetTitle("Y[cm]");
      h_EventDisplay->GetYaxis()->CenterTitle(1);
      h_EventDisplay->Draw("colztext");
      ReverseXAxis(h_EventDisplay);
      eventNumber = EVENTSPERSPILL*int(eventNumber/(EVENTSPERSPILL*nLayers)) + eventNumber%EVENTSPERSPILL;
      plotName = Form("EventDisplay_Event%d_Layer%d", eventNumber, layerNumber);
      outputPlotFullPath = (outputFolder + Form("/Detailed/")) + plotName + Form("_%06d.png", runNumber);
      plotCanvas->SaveAs(outputPlotFullPath);
    }
    delete h_EventDisplay;
  } // loop over events ends here
  
  delete inputFile;

  delete plotCanvas;
  //   }//loop over spill ends here
  // }//loop over run ends here

  // OLD:

  // char filename[50], dirname[100], histname[50];
  // int run = 2;
  // int spill = 3;

  // for(int iii=0; iii<=run; iii++){
  //   for(int jjj=1; jjj<=spill;jjj++){   
  // if(iii != 0) sprintf(filename,"test_DigiAndRechitPlotter_TB_Run%i_Spill%i.root",iii,jjj);
  // else sprintf(filename,"test_DigiAndRechitPlotter_TB_Cumulative.root");

  // const int sensorsize = 128;
  // const int Sensor_Iu = 0;
  // const int Sensor_Iv = 0;

  // HGCalTBTopology IsCellValid;
  // plotName = Form("FullLayer_RecHits_Layer%i",layerCounter);
    // // sprintf(histname,"FullLayer_RecHits_Layer%i",layerCounter);
    // TH2Poly* Overview_Rechits_Layer  = (TH2Poly*) inputFile->FindObjectAny(plotName);
    // Overview_Rechits_Layer->GetYaxis()->SetTitle("Y[cm]");
    // Overview_Rechits_Layer->GetYaxis()->CenterTitle(1);
    // Overview_Rechits_Layer->Draw("colztext");
    // ReverseXAxis(Overview_Rechits_Layer);
    // // if(iii != 0) sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_RecHits_Layer%i.png",iii,jjj,layerCounter);
    // // else sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_RecHits_Layer%i.png",layerCounter);
    // outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    // plotCanvas->SaveAs(outputPlotFullPath);
    // delete Overview_Rechits_Layer;

    // plotName = Form("FullLayer_Occupancy_Layer%i", layerCounter);
    // TH2Poly* Overview_Occupancy_Layer  = (TH2Poly*) inputFile->FindObjectAny(plotName);
    // Overview_Occupancy_Layer->GetYaxis()->SetTitle("Y[cm]");
    // Overview_Occupancy_Layer->GetYaxis()->CenterTitle(1);
    // Overview_Occupancy_Layer->Draw("colztext");
    // ReverseXAxis(Overview_Occupancy_Layer);
    // // if(iii != 0) sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_Occupancy_Layer%i.png",iii,jjj,layerCounter);
    // // else sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_Occupancy_Layer%i.png",layerCounter);
    // outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    // plotCanvas->SaveAs(outputPlotFullPath);
    // delete Overview_Occupancy_Layer;
    
    // // plotName = Form("FullLayer_RecHits_Layer%i_Summed",layerCounter);
    // // TH1F* Overview_Rechits_Layer_Summed  = (TH1F*) inputFile->FindObjectAny(plotName);
    // // Overview_Rechits_Layer_Summed->Draw();
    // // // if(iii != 0) sprintf(outputPlotFullPath,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Overview/FullLayer_RecHits_Layer%i_Summed.png",iii,jjj,layerCounter);
    // // // else sprintf(outputPlotFullPath,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Overview/FullLayer_RecHits_Layer%i_Summed.png",layerCounter);
    // // outputPlotFullPath = (outputFolder + Form("/")) + plotName + Form("_%06d.png", runNumber);
    // // plotCanvas->SaveAs(outputPlotFullPath);
    // // delete Overview_Rechits_Layer_Summed;
    
    // // // for(int iv = -7; iv < 8; iv++){
    // // //   for(int iu = -7; iu < 8; iu++){
    // // //     if(!IsCellValid.iu_iv_valid(nLayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
    // // //     plotName = Form("Cell_RecHits_u_%i_v_%i_Layer%i",iu,iv,layerCounter);
    // // //     TH1F* Cell_Hist_U_V  = (TH1F*) inputFile->FindObjectAny(plotName);
    // // //     // if(iii != 0) sprintf(outputPlotFullPath,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Run%i/Spill%i/Detailed/Cell_RecHits_u_%i_v_%i_Layer_%i.png",iii,jjj,iu,iv,layerCounter);
    // // //     // else sprintf(outputPlotFullPath,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots/Cumulative/Detailed/Cell_RecHits_u_%i_v_%i_Layer_%i.png",iu,iv,layerCounter);
    // // //     outputPlotFullPath = (outputFolder + Form("/Detailed/")) + plotName + Form("_%06d.png", runNumber);
    // // //     Cell_Hist_U_V->Draw();
    // // //     plotCanvas->SaveAs(outputPlotFullPath);
    // // //     delete Cell_Hist_U_V;
    // // //   }//loop over iu ends here
    // // // }//loop over iv ends here
}
