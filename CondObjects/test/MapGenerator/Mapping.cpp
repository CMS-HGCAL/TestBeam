#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;
ofstream fs;
bool ix_iv_valid(int ix, int iv, int sensorSize);
void swap(int (&a)[69], int (&b)[69], int (&c)[69], int size);
void Print_Map_Line(int sk, int channel, int layer, int sensor_iu, int sensor_iv, int cell_iu, int cell_iv, int type);

main(){
  const int size = 69;
  ifstream Channel_SKIROC1("Channel_Numbers_SKIROC1.txt");// Mapping of the first SKIROC(upper half in the u,v system)
  ifstream Channel_SKIROC2("Channel_Numbers_SKIROC2.txt");// Mapping of the first SKIROC(lower half in the u,v system)
  fs.open("/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/CondObjects/data/map_FNAL_1.txt");

  if(!Channel_SKIROC1) {
        cout << "Unable to open file Channel_Numbers_SKIROC1.txt";
        exit(1); // terminate with error
    }
  if(!Channel_SKIROC2) {
        cout << "Unable to open file Channel_Numbers_SKIROC2.txt";
        exit(1); // terminate with error
    }
  int sensorsize = 128;
//////Hard Coded entries for the calibration pads////////////
  int Channel_SK1_Calib = 0;
  int Channel_SK2_Calib = 0;
//Calib pad 1, closer to the central full-hex
  int IU_1_Calib = -1;
  int IV_1_Calib = 2;
//Calib pad 2, farther from the central full-hex w.r.t Calb pad 1
  int IU_2_Calib = 2;
  int IV_2_Calib = -4;
////////////////////////////////////////////////////////////

  int Channel_1[size]={0};
  int IU_1[size]={0};
  int IV_1[size]={0};
  int Channel_2[size]={0};
  int IU_2[size]={0};
  int IV_2[size]={0};
  int iterator = 1;
//Print the heading
  fs<<"# "<<"SKIROC"<<" "<<"CHANNEL"<<" | "<<"LAYER"<<" "<<"SENSOR_IX"<<" "<<"SENSOR_IV"<<" "<<"IX"<<" "<<"IV"<<" "<<"TYPE"<<endl;
//filling the El ID -> Det ID association for SKIROC 1
  for(int vvv=0;vvv<8;vvv++){
     for(int uuu=-7;uuu<8;uuu++){
        if(!ix_iv_valid(uuu,vvv,sensorsize)) continue;
        if(vvv == 0 && uuu > 0) continue;
        Channel_SKIROC1>>Channel_1[iterator] ;
        IU_1[iterator] = uuu;
        IV_1[iterator++] = vvv;
       }
    }

//filling the El ID -> Det ID association for SKIROC 1
  iterator = 1;
  for(int vvv=0;vvv> -8;vvv--){ 
     for(int uuu=-7;uuu<8;uuu++){
        if(!ix_iv_valid(uuu,vvv,sensorsize)) continue;
        if(vvv == 0 && uuu <= 0) continue;
        Channel_SKIROC2>>Channel_2[iterator];
        IU_2[iterator] = uuu;
        IV_2[iterator++] = vvv;
       }
    }

////////Sort by SKIROC channel number///////////
   swap(Channel_1,IU_1,IV_1,67);
   swap(Channel_2,IU_2,IV_2,66);
////////////////////////////////////////////////


//////////Currently we have only one sensor per layer, cell type to be implemeted////////////////////
   int sensor_iu = 0;
   int sensor_iv = 0;
   int type = 0;   
////////////////////////////////////////////////////////////////////////////////////////////////////

//Print the El ID -> Det ID map
// Due to layer rotation a cell marked u,v in odd layers goes to (u+v),-v in even layer. The logic used in this scheme is that the relative orientation of axes do not change on layer inversion in alternate layers. The ElID->DetID mapping is consistently flipped accordingly.
int hike = 0;
for(int ILayer=1;ILayer<=4;ILayer++){
   hike = 2*(ILayer - 1);
   for(int sk = 1; sk <= 2; sk++){
       if((ILayer % 2) == 1){
          if(sk == 1) Print_Map_Line(sk + hike,Channel_SK1_Calib,ILayer,sensor_iu,sensor_iv,IU_1_Calib,IV_1_Calib,type);
          if(sk == 2) Print_Map_Line(sk + hike,Channel_SK2_Calib,ILayer,sensor_iu,sensor_iv,IU_2_Calib,IV_2_Calib,type);;
       }
       else{
            if(sk == 1) Print_Map_Line(sk + hike,Channel_SK1_Calib,ILayer,sensor_iu,sensor_iv,(IU_1_Calib+IV_1_Calib),-IV_1_Calib,type);
            if(sk == 2) Print_Map_Line(sk + hike,Channel_SK2_Calib,ILayer,sensor_iu,sensor_iv,(IU_2_Calib+IV_2_Calib),-IV_2_Calib,type);
           } 
       for(int iii = 1; iii <= 67; iii++){
           if(sk == 2 && iii > 66) continue;
           if((ILayer % 2) == 1){ 
              if(sk == 1) Print_Map_Line(sk + hike,Channel_1[iii],ILayer,sensor_iu,sensor_iv,IU_1[iii],IV_1[iii],type);
              if(sk == 2) Print_Map_Line(sk + hike,Channel_2[iii],ILayer,sensor_iu,sensor_iv,IU_2[iii],IV_2[iii],type);
             }
           else{
                if(sk == 1) Print_Map_Line(sk + hike,Channel_1[iii],ILayer,sensor_iu,sensor_iv,(IU_1[iii]+IV_1[iii]),-IV_1[iii],type);
                if(sk == 2) Print_Map_Line(sk + hike,Channel_2[iii],ILayer,sensor_iu,sensor_iv,(IU_2[iii]+IV_2[iii]),-IV_2[iii],type);
               } 
           }// loop over the channels per SKIROC ends here
        }// loop over 2 SKIROCs per layer ends here
 }// loop over layer ends here


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
          }// loop over jjj ends here
       }// loop over iii ends here
   }

void Print_Map_Line(int sk, int channel, int layer, int sensor_iu, int sensor_iv, int cell_iu, int cell_iv, int type){
     fs<<"\t"<<sk<<"\t"<<channel<<"\t"<<layer<<"\t"<<sensor_iu<<"\t"<<sensor_iv<<"\t"<<cell_iu<<"\t"<<cell_iv<<"\t"<<type<<endl;
    }



