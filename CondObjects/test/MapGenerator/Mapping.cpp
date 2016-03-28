#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;
ofstream fs;
bool ix_iv_valid(int ix, int iv, int sensorSize);
void swap(int (&a)[69], int (&b)[69], int (&c)[69], int size);
void Print_Map_Line(int sk, int channel, int layer, int sensor_iu, int sensor_iv, int cell_iu, int cell_iv, int type);
int Cell_Type(int channel, int skiroc);// 0 for Full-Hex, 1 for calib pad, 2 for half-hex, 3 for mouse-bite, 4 outer calib-pad, 5 for merged-cell

//////Harcoded entries for the calibration pad////////////////////////////
const int number_of_calib_pads = 1;
int Map_Calib_Pad[number_of_calib_pads] = {0}; 
///////////////Hardcoded Map entries for mouse bites/////////////////
const int number_of_mouse_bites = 3;// Please note its 3 per sensor
int Map_Mouse_Bites_SK1[number_of_mouse_bites] = {7,38,58};
int Map_Mouse_Bites_SK2[number_of_mouse_bites] = {14,32,48};
/////////////////Hardcoded Map entries for half hex////////////////////
const int number_of_half_hex_SK1 = 5;// one less due to merges 12SR#1
const int number_of_half_hex_SK2 = 6;
int Map_Half_Hex_SK1[number_of_half_hex_SK1] = {4, 28, 45, 50, 61};
int Map_Half_Hex_SK2[number_of_half_hex_SK2] = {9, 23, 28, 38, 42, 52};
//////////////////Hardcoded entry for merged cell//////////////////////
const int number_of_merged_cells = 1;
int Map_Merged_Cells_SKI1[number_of_merged_cells] = {12};
/////////////////Hardcoded entries Outer Calib Pad////////////////////////////
const int number_of_outer_calib_pads = 1;
int Map_Merged_Outer_Calib_Pad_SK1[number_of_outer_calib_pads] = {35};
int Map_Merged_Outer_Calib_Pad_SK2[number_of_outer_calib_pads] = {29};
////////Hard Coded entries for the calibration pads////////////
//Calib pad 1, closer to the central full-hex
int IU_1_Calib = -1;
int IV_1_Calib = 2;
//Calib pad 2, farther from the central full-hex w.r.t Calb pad 1
int IU_2_Calib = 2;
int IV_2_Calib = -4;
///////////////////////////////////////////////////////////


main(){
  const int size = 69;
  ifstream Channel_SKIROC1("Channel_Numbers_SKIROC1.txt");// Mapping of the first SKIROC(upper half in the u,v system)
  ifstream Channel_SKIROC2("Channel_Numbers_SKIROC2.txt");// Mapping of the first SKIROC(lower half in the u,v system)
  fs.open("/afs/cern.ch/work/r/rchatter/CMSSW_7_6_3_patch2/src/HGCal/CondObjects/data/map_FNAL_2.txt");

  if(!Channel_SKIROC1) {
        cout << "Unable to open file Channel_Numbers_SKIROC1.txt";
        exit(1); // terminate with error
    }
  if(!Channel_SKIROC2) {
        cout << "Unable to open file Channel_Numbers_SKIROC2.txt";
        exit(1); // terminate with error
    }

  int sensorsize = 128;
  int Channel_1[size]={0};
  int IU_1[size]={0};
  int IV_1[size]={0};
  int Channel_2[size]={0};
  int IU_2[size]={0};
  int IV_2[size]={0};
  int iterator = 0;
//Print the heading
  fs<<"# "<<"SKIROC"<<" "<<"CHANNEL"<<" | "<<"LAYER"<<" "<<"SENSOR_IX"<<" "<<"SENSOR_IV"<<" "<<"IX"<<" "<<"IV"<<" "<<"TYPE"<<endl;
//filling the El ID -> Det ID association for SKIROC 1
  Channel_SKIROC1>>Channel_1[iterator] ;
  if(iterator == 0){
     IU_1[iterator] = IU_1_Calib;
     IV_1[iterator++] = IV_1_Calib;
    }   
  for(int vvv=0;vvv<8;vvv++){
     for(int uuu=-7;uuu<8;uuu++){
        if(!ix_iv_valid(uuu,vvv,sensorsize)) continue;
        if(vvv == 0 && uuu > 0) continue;
        Channel_SKIROC1>>Channel_1[iterator] ;
        IU_1[iterator] = uuu;
        IV_1[iterator++] = vvv;
       }
    }
  swap(Channel_1,IU_1,IV_1,iterator);//Sort by SKIROC channel number
//filling the El ID -> Det ID association for SKIROC 2
  iterator = 0;
  Channel_SKIROC2>>Channel_2[iterator];
  if(iterator == 0){
     IU_2[iterator] = IU_2_Calib;
     IV_2[iterator++] = IV_2_Calib;
    }
  for(int vvv=0;vvv> -8;vvv--){ 
     for(int uuu=-7;uuu<8;uuu++){
        if(!ix_iv_valid(uuu,vvv,sensorsize)) continue;
        if(vvv == 0 && uuu <= 0) continue;
        Channel_SKIROC2>>Channel_2[iterator];
        IU_2[iterator] = uuu;
        IV_2[iterator++] = vvv;
       }
    }
  swap(Channel_2,IU_2,IV_2,67);//Sort by SKIROC channel number

//////////Currently we have only one sensor per layer, cell type to be implemeted////////////////////
   int sensor_iu = 0;
   int sensor_iv = 0;
////////////////////////////////////////////////////////////////////////////////////////////////////

//Print the El ID -> Det ID map
// Due to layer rotation a cell marked u,v in odd layers goes to (u+v),-v in even layer. The logic used in this scheme is that the relative orientation of axes do not change on layer inversion in alternate layers. The ElID->DetID mapping is consistently flipped accordingly.
int hike = 0;
for(int ILayer=1;ILayer<=4;ILayer++){
   hike = 2*(ILayer - 1);
   for(int sk = 1; sk <= 2; sk++){

       for(int iii = 0; iii < 68; iii++){
           if(sk == 2 && iii > 66) continue;
           if((ILayer % 2) == 1){ 
              if(sk == 1) Print_Map_Line(sk + hike,Channel_1[iii],ILayer,sensor_iu,sensor_iv,IU_1[iii],IV_1[iii],Cell_Type(Channel_1[iii], sk));
              if(sk == 2) Print_Map_Line(sk + hike,Channel_2[iii],ILayer,sensor_iu,sensor_iv,IU_2[iii],IV_2[iii],Cell_Type(Channel_2[iii], sk));
             }
           else{
                if(sk == 1) Print_Map_Line(sk + hike,Channel_1[iii],ILayer,sensor_iu,sensor_iv,(IU_1[iii]+IV_1[iii]),-IV_1[iii],Cell_Type(Channel_1[iii], sk));
                if(sk == 2) Print_Map_Line(sk + hike,Channel_2[iii],ILayer,sensor_iu,sensor_iv,(IU_2[iii]+IV_2[iii]),-IV_2[iii],Cell_Type(Channel_1[iii], sk));
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
     for(int iii=0;iii<size-1;iii++){
        for(int jjj=iii+1; jjj<size; jjj++){
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

int Cell_Type(int channel, int skiroc){
    if((skiroc % 2) == 1){

       for(int iii=0;iii<number_of_calib_pads;iii++)
          if(channel == Map_Calib_Pad[iii]) return 1;

       for(int iii=0;iii<number_of_half_hex_SK1;iii++)
          if(channel == Map_Half_Hex_SK1[iii]) return 2;

       for(int iii=0;iii<number_of_mouse_bites;iii++) 
          if(channel == Map_Mouse_Bites_SK1[iii]) return 3;

       for(int iii=0;iii<number_of_outer_calib_pads;iii++)
          if(channel == Map_Merged_Outer_Calib_Pad_SK1[iii]) return 4;

       for(int iii=0;iii<number_of_merged_cells;iii++)
          if(channel == Map_Merged_Cells_SKI1[iii]) return 5;

       return 0;// It has to be a full hex
     }
    else{
           for(int iii=0;iii<number_of_calib_pads;iii++)
               if(channel == Map_Calib_Pad[iii]) return 1;

           for(int iii=0;iii<number_of_half_hex_SK2;iii++)
               if(channel == Map_Half_Hex_SK2[iii]) return 2;

           for(int iii=0;iii<number_of_mouse_bites;iii++) 
               if(channel == Map_Mouse_Bites_SK2[iii]) return 3;

           for(int iii=0;iii<number_of_outer_calib_pads;iii++)
               if(channel == Map_Merged_Outer_Calib_Pad_SK2[iii]) return 4;

           return 0;// It has to be a full hex
          }
  } 
