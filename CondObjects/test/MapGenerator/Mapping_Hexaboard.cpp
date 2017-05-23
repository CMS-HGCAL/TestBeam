#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <map>
using namespace std;
ofstream fs;
bool ix_iv_valid(int ix, int iv, int sensorSize);
void swap(int (&a)[69], int (&b)[69], int (&c)[69], int size);
void Print_Map_Line(int sk, int channel, int layer, int sensor_iu, int sensor_iv, int cell_iu, int cell_iv, int type);
int Cell_Type(int channel, int skiroc);// 0 for Full-Hex, 1 for calib pad, 2 for half-hex, 3 for mouse-bite, 4 outer calib-pad, 5 for merged-cell

//////Harcoded entries for the calibration pad////////////////////////////
/*
const int number_of_calib_pads = 1;
int Map_Calib_Pad_SK1[number_of_calib_pads] = {22};
int Map_Calib_Pad_SK2[number_of_calib_pads] = {22}; 
*/
///////////////Hardcoded Map entries for mouse bites/////////////////////
const int number_of_mouse_bites_SK1 = 1;
const int number_of_mouse_bites_SK2 = 2;
const int number_of_mouse_bites_SK3 = 1;
const int number_of_mouse_bites_SK4 = 2;
int Map_Mouse_Bites_SK1[number_of_mouse_bites_SK1] = {46};
int Map_Mouse_Bites_SK2[number_of_mouse_bites_SK2] = {20, 60};
int Map_Mouse_Bites_SK3[number_of_mouse_bites_SK3] = {12};
int Map_Mouse_Bites_SK4[number_of_mouse_bites_SK4] = {34, 16};
/////////////////Hardcoded Map entries for half hex////////////////////
const int number_of_half_hex = 3;
int Map_Half_Hex_SK1[number_of_half_hex] = {10, 56, 50};
int Map_Half_Hex_SK2[number_of_half_hex] = {12, 10, 48};
int Map_Half_Hex_SK3[number_of_half_hex] = {16, 62, 52};
int Map_Half_Hex_SK4[number_of_half_hex] = {26, 22, 60};
//////////////////Hardcoded entry for merged cell//////////////////////
/*
const int number_of_merged_cells = 2;
int Map_Merged_Cells_SKI1[number_of_merged_cells] = {0,60};
const int number_of_merged_cells_SKI2 = 2;
int Map_Merged_Cells_SKI2[number_of_merged_cells] = {0,60};
*/
/////////////////Hardcoded entries Outer Calib Pad////////////////////////////
const int number_of_outer_calib_pads = 1;
int Map_Merged_Outer_Calib_Pad_SK1[number_of_outer_calib_pads] = {2};
int Map_Merged_Outer_Calib_Pad_SK3[number_of_outer_calib_pads] = {26};
///////////////////////////////////////////////////////////


main(){


/*
  ifstream MapFile("Mapping.txt");
  int Cell, Channel;
  map<int,int> SkirocToCell;
  map<int,int>::iterator Ski_It;
  int MapFileCounter=1;
  while(!MapFile.eof() && MapFileCounter<=127){
      MapFile>>Cell;
      MapFile>>Channel;
      SkirocToCell[Cell] = Channel;
     }
*/

  const int size = 69;
  ifstream Channel_SKIROC1("Channel_Numbers_SKIROC1_Hexaboard_Updated.txt");// Mapping of the first SKIROC(upper half in the u,v system)
  ifstream Channel_SKIROC2("Channel_Numbers_SKIROC2_Hexaboard_Updated.txt");// Mapping of the first SKIROC(lower half in the u,v system)
  ifstream Channel_SKIROC3("Channel_Numbers_SKIROC3_Hexaboard_Updated.txt");// Mapping of the first SKIROC(lower half in the u,v system)
  ifstream Channel_SKIROC4("Channel_Numbers_SKIROC4_Hexaboard_Updated.txt");// Mapping of the first SKIROC(lower half in the u,v system)
//  fs.open("/afs/cern.ch/work/r/rchatter/May_Test_Beam_Test/CMSSW_8_0_1/src/HGCal/CondObjects/data/map_FNAL_SB1.txt");
  fs.open("map_CERN_Hexaboard_28Layers.txt");

  if(!Channel_SKIROC1) {
        cout << "Unable to open file Channel_Numbers_SKIROC1.txt";
        exit(1); // terminate with error
    }
  if(!Channel_SKIROC2) {
        cout << "Unable to open file Channel_Numbers_SKIROC2.txt";
        exit(1); // terminate with error
    }
  if(!Channel_SKIROC3) {
        cout << "Unable to open file Channel_Numbers_SKIROC3.txt";
        exit(1); // terminate with error
    }
  if(!Channel_SKIROC4) {
        cout << "Unable to open file Channel_Numbers_SKIROC4.txt";
        exit(1); // terminate with error
    }

  int sensorsize = 128;
  int Channel_1[size]={0};
  int IU_1[size]={0};
  int IV_1[size]={0};
  int Channel_2[size]={0};
  int IU_2[size]={0};
  int IV_2[size]={0};
  int Channel_3[size]={0};
  int IU_3[size]={0};
  int IV_3[size]={0};
  int Channel_4[size]={0};
  int IU_4[size]={0};
  int IV_4[size]={0};
  int iterator = 0;
  int SKI_Channels[4] = {0};

//Print the heading
  fs<<"# "<<"SKIROC"<<" "<<"CHANNEL"<<" | "<<"LAYER"<<" "<<"SENSOR_IX"<<" "<<"SENSOR_IV"<<" "<<"IX"<<" "<<"IV"<<" "<<"TYPE"<<endl;
  for(int vvv = -7; vvv <= 7; vvv++){
     for(int uuu = -7; uuu < 8; uuu++){
        if(!ix_iv_valid(uuu, vvv, sensorsize)) continue;

        if(vvv == -7 && uuu !=  4) continue;
        if(vvv == -6 && uuu >  5) continue;
        if(vvv == -5 && uuu >  4) continue;
        if(vvv == -4 && (uuu == -3 || uuu >  3)) continue;
        if(vvv == -3 && (uuu == -4 || uuu >  3)) continue;
        if(vvv == -2 && (uuu < -2 || uuu >  2)) continue;
        if(vvv == -1 && (uuu < -1 || uuu >  0)) continue;
	if(vvv >= 0) continue;
        Channel_SKIROC1>>Channel_1[iterator];
        IU_1[iterator] = uuu;
        IV_1[iterator++] = vvv;
       }
    }

  swap(Channel_1,IU_1,IV_1,iterator);
  SKI_Channels[0] = iterator;
	
  iterator = 0;

  for(int vvv = -7; vvv < 8; vvv++){ 
     for(int uuu = -7; uuu < 8; uuu++){
        if(!ix_iv_valid(uuu, vvv, sensorsize)) continue;
	if(vvv < -3) continue;
        if(vvv ==  -3 && uuu != -4 ) continue;
	if(vvv ==  -2 && uuu > -3 ) continue;
	if(vvv ==  -1 && uuu > -2 ) continue;
	if(vvv ==  0 && uuu > -1 ) continue;
        if(vvv ==  1 && uuu > -1) continue;
        if(vvv ==  2 && uuu > -3) continue;
        if(vvv ==  3 && (uuu == -7 || uuu > -3)) continue;
        if(vvv ==  4 && uuu > -4) continue;
        if(vvv ==  5 && uuu > -5) continue;
        if(vvv >  5) continue;
        Channel_SKIROC2>>Channel_2[iterator];
        IU_2[iterator] = uuu;
        IV_2[iterator++] = vvv;
       }
    }
  swap(Channel_2,IU_2,IV_2,iterator);
  SKI_Channels[1] = iterator;

  iterator = 0;

  for(int vvv = -7; vvv < 8; vvv++){
     for(int uuu = -7; uuu < 8; uuu++){
        if(!ix_iv_valid(uuu, vvv, sensorsize)) continue;
	if(vvv < 0) continue;
        if(vvv ==  0 &&  uuu != 0 ) continue;
        if(vvv ==  1 && (uuu < 0 || uuu >  1)) continue;
        if(vvv ==  2 && (uuu < -2 || uuu >  2)) continue;
        if(vvv ==  3 && (uuu < -2 || uuu >  2)) continue;
        if(vvv ==  4 && (uuu < -3 || uuu >  2)) continue;
        if(vvv ==  5 && (uuu < -4 || uuu >  1)) continue;
        if(vvv ==  6 && uuu < -5) continue;
        if(vvv ==  7 && uuu != -3) continue;
        Channel_SKIROC3>>Channel_3[iterator];
        IU_3[iterator] = uuu;
        IV_3[iterator++] = vvv;
       }
    }
  swap(Channel_3,IU_3,IV_3,iterator);
  SKI_Channels[2] = iterator;

  iterator = 0;

  for(int vvv = -7; vvv <= 8; vvv++){
     for(int uuu = -7; uuu < 8; uuu++){
        if(!ix_iv_valid(uuu, vvv, sensorsize)) continue;
	if(vvv < -5) continue;
        if(vvv == -5 && uuu <  5) continue;
        if(vvv == -4 && uuu <  4) continue;
        if(vvv == -3 && (uuu == 7 || uuu <  4)) continue;
        if(vvv == -2 && uuu <  3) continue;
        if(vvv == -1 && uuu <  1) continue;
        if(vvv ==  0 && uuu <  1) continue;
	if(vvv ==  1 && uuu < 2) continue;
	if(vvv ==  2 && uuu < 3) continue;
	if(vvv ==  3 && uuu != 3) continue;
	if(vvv ==  4 && uuu != 3) continue;
	if(vvv > 4) continue;
        Channel_SKIROC4>>Channel_4[iterator];
        IU_4[iterator] = uuu;
        IV_4[iterator++] = vvv;
       }
    }
  swap(Channel_4,IU_4,IV_4,iterator);
  SKI_Channels[3] = iterator;
//////////Currently we have only one sensor per layer, cell type to be implemeted////////////////////
   int sensor_iu = 0;
   int sensor_iv = 0;
////////////////////////////////////////////////////////////////////////////////////////////////////
//Print the El ID -> Det ID map
// Due to layer rotation a cell marked u,v in odd layers goes to (u+v),-v in even layer. The logic used in this scheme is that the relative orientation of axes do not change on layer inversion in alternate layers. The ElID->DetID mapping is consistently flipped accordingly.
int hike = 0;
for(int ILayer = 1; ILayer <= 28; ILayer++){
   hike = 4*(ILayer - 1);
   for(int sk = 1; sk <= 4; sk++){

       for(int iii = 0; iii < 68; iii++){
           if( iii >= SKI_Channels[sk - 1] ) continue;
 if(ILayer%2 == 0){
                   if(sk == 1) Print_Map_Line(sk + hike,Channel_1[iii],ILayer,sensor_iu,sensor_iv,IU_1[iii],IV_1[iii],Cell_Type(Channel_1[iii], sk));
                   if(sk == 2) Print_Map_Line(sk + hike,Channel_2[iii],ILayer,sensor_iu,sensor_iv,IU_2[iii],IV_2[iii],Cell_Type(Channel_2[iii], sk));
		   if(sk == 3) Print_Map_Line(sk + hike,Channel_3[iii],ILayer,sensor_iu,sensor_iv,IU_3[iii],IV_3[iii],Cell_Type(Channel_3[iii], sk));
		   if(sk == 4) Print_Map_Line(sk + hike,Channel_4[iii],ILayer,sensor_iu,sensor_iv,IU_4[iii],IV_4[iii],Cell_Type(Channel_4[iii], sk));
                   }
              else{
                    if(sk == 1) Print_Map_Line(sk + hike,Channel_1[iii],ILayer,sensor_iu,sensor_iv,(IU_1[iii]+IV_1[iii]),-IV_1[iii],Cell_Type(Channel_1[iii], sk));
                    if(sk == 2){
				 Print_Map_Line(sk + hike,Channel_2[iii],ILayer,sensor_iu,sensor_iv,(IU_2[iii]+IV_2[iii]),-IV_2[iii],Cell_Type(Channel_2[iii], sk));
			}
		    if(sk == 3) Print_Map_Line(sk + hike,Channel_3[iii],ILayer,sensor_iu,sensor_iv,(IU_3[iii]+IV_3[iii]),-IV_3[iii],Cell_Type(Channel_3[iii], sk));
		    if(sk == 4) Print_Map_Line(sk + hike,Channel_4[iii],ILayer,sensor_iu,sensor_iv,(IU_4[iii]+IV_4[iii]),-IV_4[iii],Cell_Type(Channel_4[iii], sk));
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
    if(skiroc  == 1){
       for(int iii=0;iii<number_of_half_hex;iii++)
          if(channel == Map_Half_Hex_SK1[iii]) return 2;

       for(int iii=0;iii<number_of_mouse_bites_SK1;iii++) 
          if(channel == Map_Mouse_Bites_SK1[iii]) return 3;
       return 0;// It has to be a full hex
     }
    else if(skiroc  == 2){
           for(int iii=0;iii<number_of_half_hex;iii++)
               if(channel == Map_Half_Hex_SK2[iii]) return 2;

           for(int iii=0;iii<number_of_mouse_bites_SK2;iii++) 
               if(channel == Map_Mouse_Bites_SK2[iii]) return 3;

           for(int iii=0;iii<number_of_outer_calib_pads;iii++)
               if(channel == Map_Merged_Outer_Calib_Pad_SK2[iii]) return 4;
           return 0;
          }
   else if(skiroc == 3){
           for(int iii=0;iii<number_of_half_hex;iii++)
               if(channel == Map_Half_Hex_SK3[iii]) return 2;

           for(int iii=0;iii<number_of_mouse_bites_SK3;iii++)
               if(channel == Map_Mouse_Bites_SK3[iii]) return 3;
           return 0;
   	}
  else{
           for(int iii=0;iii<number_of_half_hex;iii++) 
               if(channel == Map_Half_Hex_SK4[iii]) return 2;

           for(int iii=0;iii<number_of_mouse_bites_SK4;iii++)
               if(channel == Map_Mouse_Bites_SK4[iii]) return 3;

           for(int iii=0;iii<number_of_outer_calib_pads;iii++)
               if(channel == Map_Merged_Outer_Calib_Pad_SK4[iii]) return 4;
           return 0;
	} 
  } 
