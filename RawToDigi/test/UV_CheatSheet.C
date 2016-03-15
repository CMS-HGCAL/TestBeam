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
   newaxis->Draw();
}


void UV_CheatSheet(){
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TFile *F = TFile::Open("test_ElectronicsMap_TB.root");
F->cd("plot");
TH2Poly* h = (TH2Poly*) F->FindObjectAny("FullLayer_Layer1_U_V");
h->Draw("colztext");
ReverseXAxis(h);

HGCalTBTopology IsCellValid;
HGCalTBCellVertices TheCell;
std::pair<double, double> CellCentreXY;
int sensorsize = 128;
int Sensor_Iu = 0;
int Sensor_Iv = 0;
int nlayers = 1;

TLatex *t = new TLatex();
t->SetTextFont(32);
t->SetTextColor(1);
t->SetTextSize(0.03);
t->SetTextAlign(12);
char U_V[50];


    for(int iv = -7; iv < 8; iv++){
        for(int iu = -7; iu < 8; iu++){
            if(!IsCellValid.iu_iv_valid(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize)) continue;
            CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots(nlayers, Sensor_Iu, Sensor_Iv, iu, iv, sensorsize);
            double iux = (CellCentreXY.first < 0 ) ? (CellCentreXY.first + 0.0001) : (CellCentreXY.first - 0.0001) ;
            double iyy = (CellCentreXY.second < 0 ) ? (CellCentreXY.second + 0.0001) : (CellCentreXY.second - 0.0001);
            sprintf(U_V,"%i,%i",iu,iv);
            t->DrawLatex(iux,0.975*iyy,U_V);
           }
       }

}


