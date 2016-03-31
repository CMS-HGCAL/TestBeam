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
int adcs= 2;
const int layers = 4;
const int sensorsize = 128;
const int Sensor_Iu = 0;
const int Sensor_Iv = 0;

HGCalTBTopology IsCellValid;
TFile *F;

for(int iii=1; iii<=run; iii++){
	for(int kkk=1; kkk<= layers; kkk++){   
		
		for(int jjj=0;jjj<adcs;jjj++){
			sprintf(filename,"test_DigiAndRechitPlotter_TB_Run%i.root",iii);
			F = TFile::Open(filename);
			F->cd("hgcaltbdigisplotter");

			sprintf(histname,"FullLayer_ADC%i_Layer%i",jjj,kkk);
			TH2Poly* Overview_Digis_Layer  = (TH2Poly*) F->FindObjectAny(histname);
                        Overview_Digis_Layer->Scale(1./512,"");
			Overview_Digis_Layer->GetYaxis()->SetTitle("Y[cm]");
			Overview_Digis_Layer->GetYaxis()->CenterTitle(1);
                        Overview_Digis_Layer->Draw("colztext");
			ReverseXAxis(Overview_Digis_Layer);
                	sprintf(dirname,"/home/rchatter/shervinTest/CMSSW_7_6_3_patch2/src/HGCal/Digis/Run%i/FullLayer_Digis_ADC%i_Layer%i.png",iii,jjj,kkk);
                        c1->SaveAs(dirname);
 
                        sprintf(histname,"FullLayer_ADC%i_Layer%i_summed",jjj,kkk);
			TH1F* Overview_Digis_Layer_Summed  = (TH1F*) F->FindObjectAny(histname);
                        Overview_Digis_Layer_Summed->Draw();
                        sprintf(dirname,"/home/rchatter/shervinTest/CMSSW_7_6_3_patch2/src/HGCal/Digis/Run%i/FullLayer_Digis_ADC%i_Layer%i_Summed.png",iii,jjj,kkk);
                        c1->SaveAs(dirname);
                        sprintf(histname,"FullLayer_ADC%i_Layer%i_profile",jjj,kkk);
                        TProfile* Overview_Digis_Layer_Profile  = (TProfile*) F->FindObjectAny(histname);
                        Overview_Digis_Layer_Profile->Draw();
                        sprintf(dirname,"/home/rchatter/shervinTest/CMSSW_7_6_3_patch2/src/HGCal/Digis/Run%i/FullLayer_Digis_ADC%i_Layer%i_Profile.png",iii,jjj,kkk);
                        c1->SaveAs(dirname);  
	            }// loop over adcs ends here

		for(int iv = -7; iv < 8; iv++){
			for(int iu = -7; iu < 8; iu++){
				for(int jjj=0;jjj<adcs;jjj++){
					if(!IsCellValid.iu_iv_valid(layers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
					sprintf(histname,"Cell_u_%i_v_%i_ADC%i_Layer%i",iu,iv,jjj,kkk);
					TH1F* Cell_Hist_U_V  = (TH1F*) F->FindObjectAny(histname);
					sprintf(dirname,"/home/rchatter/shervinTest/CMSSW_7_6_3_patch2/src/HGCal/Digis/Run%i/Detailed/Cell_u_%i_v_%i_ADC%i_Layer_%i.png",iii,iu,iv,jjj,kkk);
					Cell_Hist_U_V->Draw();
					c1->SaveAs(dirname);
                                    } 
                            }//loop over iu ends here
                   }//loop over iv ends here

	}//loop over layer ends here
   }//loop over run ends her

}


