#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

bool ix_iv_valid(int ix, int iv, int sensorSize);
void swap(int (&a)[69], int (&b)[69], int (&c)[69], int size);


main(){

  const int size = 69;
  ofstream fs;
  ifstream Channel_SKIROC1("Channel_Numbers_SKIROC1.txt");// Mapping of the first SKIROC(upper half in the u,v system)
  ifstream Channel_SKIROC2("Channel_Numbers_SKIROC2.txt");// Mapping of the first SKIROC(lower half in the u,v system)
  if(!Channel_SKIROC1) {
        cout << "Unable to open file Channel_Numbers_SKIROC1.txt";
        exit(1); // terminate with error
    }
  if(!Channel_SKIROC2) {
        cout << "Unable to open file Channel_Numbers_SKIROC2.txt";
        exit(1); // terminate with error
    }
  int sensorsize = 128;
  int Layer = 1;
  int SKIROC = 1;
//////Hard Coded entries for the calibration pads////////////
  int Channel_SK1_Calib = 0;
  int Channel_SK2_Calib = 0;
//Calib pad 1, closer to the central full-hex
  int IX_1_Calib = -1;
  int IV_1_Calib = 2;
//Calib pad 2, farther from the central full-hex w.r.t Calb pad 1
  int IX_2_Calib = 2;
  int IV_2_Calib = -4;
////////////////////////////////////////////////////////////

  int Channel_1[69]={0};
  int IX_1[69]={0};
  int IV_1[69]={0};
  int Channel_2[69]={0};
  int IX_2[69]={0};
  int IV_2[69]={0};
  int iterator = 1;
  fs.open("/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/CondObjects/data/map_FNAL.txt");
  fs<<"# "<<"SKIROC"<<" "<<"CHANNEL"<<" | "<<"LAYER"<<" "<<"SENSOR_IX"<<" "<<"SENSOR_IV"<<" "<<"IX"<<" "<<"IV"<<" "<<"TYPE"<<endl;
//  fs<<"Layer"<<"\t"<<"SKIROC"<<"\t"<<"Channel"<<"\t"<<"ix"<<"\t"<<"iv"<<endl;
  SKIROC = 1;
  for(int vvv=0;vvv<8;vvv++){
     for(int iii=-7;iii<8;iii++){
        if(!ix_iv_valid(iii,vvv,sensorsize)) continue;
        if(vvv == 0 && iii > 0) continue;
        Channel_SKIROC1>>Channel_1[iterator] ;
        IX_1[iterator] = iii;
        IV_1[iterator++] = vvv;
       }
    }

cout<<endl<<" Iterator1 = "<<iterator<<endl;

  SKIROC = 2;
  iterator = 1;
  for(int vvv=0;vvv> -8;vvv--){ 
     for(int iii=-7;iii<8;iii++){
        if(!ix_iv_valid(iii,vvv,sensorsize)) continue;
        if(vvv == 0 && iii <= 0) continue;
        Channel_SKIROC2>>Channel_2[iterator];
        IX_2[iterator] = iii;
        IV_2[iterator++] = vvv;
       }
    }


   swap(Channel_1,IX_1,IV_1,67);
   swap(Channel_2,IX_2,IV_2,66);
   SKIROC = 1;
   int sensor_ix = 0;
   int sensor_iv = 0;
   int type = 0;   


  for(int sk =1; sk<=2;sk++){
     if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_SK1_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_1_Calib<<"\t"<<IV_1_Calib<<"\t"<<type<<endl;
     if(sk == 2) fs<<"\t"<<sk<<"\t"<<Channel_SK2_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_2_Calib<<"\t"<<IV_2_Calib<<"\t"<<type<<endl;
    for(int iii = 1; iii<=67; iii++){
     if(sk == 2 && iii > 66) continue;
     if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_1[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_1[iii]<<"\t"<<IV_1[iii]<<"\t"<<type<<endl;
     else fs<<"\t"<<sk<<"\t"<<Channel_2[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_2[iii]<<"\t"<<IV_2[iii]<<"\t"<<type<<endl;
       }
   }

Layer = 2;
  for(int sk =1; sk<=2;sk++){
     if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_SK1_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_1_Calib<<"\t"<<IV_1_Calib<<"\t"<<type<<endl;
     if(sk == 2) fs<<"\t"<<sk<<"\t"<<Channel_SK2_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_2_Calib<<"\t"<<IV_2_Calib<<"\t"<<type<<endl;
    for(int iii = 1; iii<=67; iii++){
     if(sk == 2 && iii > 66) continue;
     if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_1[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_1[iii]<<"\t"<<IV_1[iii]<<"\t"<<type<<endl;
     else fs<<"\t"<<sk<<"\t"<<Channel_2[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_2[iii]<<"\t"<<IV_2[iii]<<"\t"<<type<<endl;
       }
   }

Layer = 3;
  for(int sk =1; sk<=2;sk++){
     if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_SK1_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_1_Calib<<"\t"<<IV_1_Calib<<"\t"<<type<<endl;
     if(sk == 2) fs<<"\t"<<sk<<"\t"<<Channel_SK2_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_2_Calib<<"\t"<<IV_2_Calib<<"\t"<<type<<endl;
    for(int iii = 1; iii<=67; iii++){
     if(sk == 2 && iii > 66) continue;
     if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_1[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_1[iii]<<"\t"<<IV_1[iii]<<"\t"<<type<<endl;
     else fs<<"\t"<<sk<<"\t"<<Channel_2[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<IX_2[iii]<<"\t"<<IV_2[iii]<<"\t"<<type<<endl;
       }
   }

Layer = 4;
  for(int sk =1; sk<=2;sk++){
    if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_SK1_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_1_Calib<<"\t"<<IV_1_Calib<<"\t"<<type<<endl;
     if(sk == 2) fs<<"\t"<<sk<<"\t"<<Channel_SK2_Calib<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_2_Calib<<"\t"<<IV_2_Calib<<"\t"<<type<<endl;
    for(int iii = 1; iii<=67; iii++){
     if(sk == 2 && iii > 66) continue;
     if(sk == 1) fs<<"\t"<<sk<<"\t"<<Channel_1[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_1[iii]<<"\t"<<IV_1[iii]<<"\t"<<type<<endl;
     else fs<<"\t"<<sk<<"\t"<<Channel_2[iii]<<"\t"<<Layer<<"\t"<<sensor_ix<<"\t"<<sensor_iv<<"\t"<<-IX_2[iii]<<"\t"<<IV_2[iii]<<"\t"<<type<<endl;
       }
   }



}//main ends here


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

void swap(int (&a)[69], int (&b)[69], int (&c)[69], int size){
     int tmp1 = 0, tmp2=0, tmp3=0;
     for(int iii=1;iii<=size;iii++){
        for(int jjj=iii+1; jjj<=size; jjj++){
           if(a[iii] > a[jjj]){
              tmp1 = a[iii]; tmp2 = b[iii]; tmp3 = c[iii];
              a[iii] = a[jjj]; b[iii] = b[jjj]; c[iii] = c[jjj]; 
              a[jjj] = tmp1; b[jjj] = tmp2; c[jjj] = tmp3;
            }
          }
       }
   }
