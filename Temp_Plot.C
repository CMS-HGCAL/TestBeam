void Temp_Plot(){
gStyle->SetOptStat(0);
gStyle->SetPaintTextFormat("3.1f");
gStyle->SetOptTitle(0);
char histname[200],dirname[200];

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 900);
  c1->Range(0, 0, 1, 1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);

  TFile* F = TFile::Open("test_DigiAndRechitPlotter_TB_8210_Ped.root");
  F->cd("hgcaltbdigisplotter_new");

  sprintf(histname,"FullLayer_ADC1_Layer1_Event75");
  TH2Poly* Overview_Digis_Layer  = (TH2Poly*) F->FindObjectAny(histname);
  Overview_Digis_Layer->SetMinimum(-100.);
  Overview_Digis_Layer->GetYaxis()->SetTitle("Y[cm]");
  Overview_Digis_Layer->GetYaxis()->CenterTitle(1);
  Overview_Digis_Layer->GetXaxis()->SetTitle("X[cm]");
  Overview_Digis_Layer->GetXaxis()->CenterTitle(1);

  Overview_Digis_Layer->Fill(-1,5,3.038);
  Overview_Digis_Layer->Fill(-3.5,5.5,-3.0);
  Overview_Digis_Layer->Fill(3.5,5.5,3.3);
  Overview_Digis_Layer->Fill(-6.5,-0.5,0.12);
  Overview_Digis_Layer->Fill(-6.5,0.5,0.12);
  Overview_Digis_Layer->Fill(6.5,-0.5,-3.66);
  Overview_Digis_Layer->Fill(-3.75,-5.5,4.0);
  Overview_Digis_Layer->Fill(3,-5.75,8.3);



  Overview_Digis_Layer->Draw("colztext");
  ReverseXAxis(Overview_Digis_Layer);
  sprintf(dirname,"/home/rchatter/shervinTest/CMSSW_7_6_3_patch2/src/HGCal/FullLayer_Digis_ADC0_Layer1_Run8210_TableDown_PedSubBySki_Event75.png");
  c1->SaveAs(dirname);


}
