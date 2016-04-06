#include <fstream>
#include <stdlib.h>
using namespace std;

bool ix_iv_valid(int ix, int iv, int sensorSize);

void Pedestal_Writer_Paul(){
     TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
     TFile* F = new TFile("test_Ped_8272.root");
     F->cd("hgcaltbdigisplotter");
     double channels[128];   
     double Mean[128];
     double RMS[128];
     char name[100];
     int nsample = 1;
     int Code = 0;
     int Type = 0;
     int Layer = 1;
     int SENSOR_IX = 0;
     int SENSOR_IV = 0;
     int sensorsize = 128;
     TH1F* h[1000];
     int counter = 0;
     ofstream fs;
     fs.open("/home/rchatter/shervinTest/CMSSW_7_6_3_patch2/src/HGCal/Ped_HighGain_Ski_Channel_8272.txt");
     fs<<"SKI CHANNEL Ped_Mean Ped_RMS"<<endl;

     int Counter = 0;
     for(int iii = 1; iii <= 2; iii++) {
         for(int chan = 0; chan < 64; chan++) {
             sprintf(name, "Ski_%i_Channel_%i_ADC%i",iii,chan,nsample);
             cout<<endl<<name<<endl;
             h[counter] = (TH1F*) F->FindObjectAny(name);
             fs<<iii<<" "<<chan<<" "<<h[counter]->GetMean()<<" "<<h[counter]->GetRMS()<<endl;
               }
//             cout<<endl<<Code<<"\t"<<Layer<<"\t"<<SENSOR_IX<<"\t"<<SENSOR_IV<<"\t"<<iu<<"\t"<<iv<<"\t"<<h[counter++]->GetMean()<<endl;
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

