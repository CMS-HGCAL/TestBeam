#include <fstream>
#include <stdlib.h>
using namespace std;

bool ix_iv_valid(int ix, int iv, int sensorSize);

void Pedestal_Writer(){

     TFile* F = new TFile("test_DigiAndRechitPlotter_TB_8272.root");
     F->cd("hgcaltbdigisplotter");
     char name[100];
     int nsample =1;
     int Code = 0;
     int Type = 0;
     int Layer = 1;
     int SENSOR_IX = 0;
     int SENSOR_IV = 0;
     int sensorsize = 128;
     TH1F* h[1000];
     int counter = 0;
     ofstream fs;
     fs.open("/home/rchatter/shervinTest/CMSSW_7_6_3_patch2/src/HGCal/Ped_HighGain_Test_8272.txt");
     fs<<"SCHEME_CODE 0"<<endl;
     fs<<"# CODE  LAYER SENSOR_IX SENSOR_IV  IX  IV TYPE  VALUE"<<endl;


     for(int iv = -7; iv < 8; iv++) {
         for(int iu = -7; iu < 8; iu++) {
             if(!ix_iv_valid(iu,iv,sensorsize)) continue;
             sprintf(name, "Cell_u_%i_v_%i_ADC%i_Layer%i", iu, iv, nsample, Layer);
             cout<<endl<<name<<endl;
             h[counter] = (TH1F*) F->FindObjectAny(name);
             fs<<" "<<Code<<" "<<Layer<<" "<<SENSOR_IX<<" "<<SENSOR_IV<<" "<<iu<<" "<<iv<<" "<<" "<<Type<<" "<<h[counter]->GetMean()<<endl;
//             cout<<endl<<Code<<"\t"<<Layer<<"\t"<<SENSOR_IX<<"\t"<<SENSOR_IV<<"\t"<<iu<<"\t"<<iv<<"\t"<<h[counter++]->GetMean()<<endl;
            }
         } 
     fs.close();
}


bool ix_iv_valid(int ix, int iv, int sensorSize){
  int aiv=abs(iv);
  int ixc=(iv<0)?(-ix):(ix);
  if(sensorSize==128) {
     if (iv==0) return (ix>=-5 && ix<=5);
     else if (aiv==1) return (ixc>=-6 && ixc<=5);
     else if (aiv==2) return (ixc>=-6 && ixc<=4);
     else if (aiv==3) return (ixc>=-7 && ixc<=4);
     else if (aiv==4) return (ixc>=-7 && ixc<=3);
     else if (aiv==5) return (ixc>=-6 && ixc<=1);
     else if (aiv==6) return (ixc>=-5 && ixc<=-1);
     else if (aiv==7) return (ixc==-3 || ixc==-4);
     else return false;
    }
   else return false;
}

