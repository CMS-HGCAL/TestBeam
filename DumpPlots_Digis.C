#include <memory>
#include <iostream>
#include "TH2Poly.h"
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


void DumpPlots_Digis(){
gStyle->SetOptStat(0);
//gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->Range(0, 0, 1, 1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);


char filename[200], dirname[200], histname[200];
int run = 1;
int spill = 1;
const int layers = 4;
const int sensorsize = 128;
const int Sensor_Iu = 0;
const int Sensor_Iv = 0;

HGCalTBTopology IsCellValid;
TFile *F;
TH2Poly *Overview_Digis_Layer_Cum;
TH1F *Overview_Digis_Layer_Summed_Cum, *Cell_Hist_U_V_Cum;

for(int iii=1; iii<=run; iii++){
	for(int kkk=1; kkk<= layers; kkk++){   
		
		for(int jjj=1;jjj<=spill;jjj++){
			sprintf(filename,"test_DigiAndRechitPlotter_TB_Run%i_Spill%i.root",iii,jjj);
			F = TFile::Open(filename);
			F->cd("hgcaltbdigisplotter");

			sprintf(histname,"FullLayer_Sample0_Layer%i",kkk);
			TH2Poly* Overview_Digis_Layer  = (TH2Poly*) F->FindObjectAny(histname);
			Overview_Digis_Layer->GetYaxis()->SetTitle("Y[cm]");
			Overview_Digis_Layer->GetYaxis()->CenterTitle(1);
			if(jjj == 1) Overview_Digis_Layer_Cum = (TH2Poly*) Overview_Digis_Layer->Clone("Overview_Digis_Layer_Cum");
//			else Overview_Digis_Layer_Cum->Add(Overview_Digis_Layer);
                        Overview_Digis_Layer->Draw("colztext");
			ReverseXAxis(Overview_Digis_Layer);
                	sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots_WebPage/Digis/Run%i/Spill%i/FullLayer_Digis_Layer%i.png",iii,jjj,kkk);
                        c1->SaveAs(dirname);
 
                        sprintf(histname,"FullLayer_Sample0_Layer%i_summed",kkk);
			TH1F* Overview_Digis_Layer_Summed  = (TH1F*) F->FindObjectAny(histname)->Clone("Overview_Digis_Layer_Summed");
			if(jjj == 1) Overview_Digis_Layer_Summed_Cum = (TH1F*) Overview_Digis_Layer_Summed->Clone("Overview_Digis_Layer_Summed_Cum");
			else Overview_Digis_Layer_Summed_Cum->Add(Overview_Digis_Layer_Summed);
                        Overview_Digis_Layer_Summed->Draw();
                        sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots_WebPage/Digis/Run%i/Spill%i/FullLayer_Digis_Layer%i_Summed.png",iii,jjj,kkk);
                        c1->SaveAs(dirname);
			delete Overview_Digis_Layer,Overview_Digis_Layer_Summed;
	            }// loop over spill ends here
		sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots_WebPage/Digis/Run%i/Cumulative/FullLayer_Digis_Layer%i.png",iii,kkk);
		Overview_Digis_Layer_Cum->Draw("colztext");
		ReverseXAxis(Overview_Digis_Layer_Cum);
		c1->SaveAs(dirname);
		sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots_WebPage/Digis/Run%i/Cumulative/FullLayer_Digis_Layer%i_Summed.png",iii,kkk);
		Overview_Digis_Layer_Summed_Cum->Draw();
		c1->SaveAs(dirname);

		for(int iv = -7; iv < 8; iv++){
			for(int iu = -7; iu < 8; iu++){
				for(int jjj=1;jjj<=spill;jjj++){
					if(!IsCellValid.iu_iv_valid(layers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
					sprintf(histname,"Cell_u_%i_v_%i_Sample0_Layer%i",iu,iv,kkk);
					TH1F* Cell_Hist_U_V  = (TH1F*) F->FindObjectAny(histname);
					if(jjj == 1) Cell_Hist_U_V_Cum = (TH1F*) Cell_Hist_U_V->Clone("Cell_Hist_U_V_Cum");
					else Cell_Hist_U_V_Cum->Add(Cell_Hist_U_V);
					sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots_WebPage/Digis/Run%i/Spill%i/Detailed/Cell_u_%i_v_%i_Layer_%i.png",iii,jjj,iu,iv,kkk);
					Cell_Hist_U_V->Draw();
					c1->SaveAs(dirname);
                                    } 
//				sprintf(dirname,"/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/DQM_Plots_WebPage/Digis/Run%i/Cumulative/Detailed/Cell_u_%i_v_%i_Layer_%i.png",iii,iu,iv,kkk);
//					Cell_Hist_U_V_Cum->Draw();
//					c1->SaveAs(dirname);
                            }//loop over iu ends here
                   }//loop over iv ends here

	}//loop over layer ends here
   }//loop over run ends her


}


