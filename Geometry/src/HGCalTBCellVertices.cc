#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include <stdlib.h>
#include "math.h"
#include <iostream>
#include <iomanip>
#define PI 3.14159265
using namespace std;

HGCalTBTopology Top;
double delta = 0.00001;//needed for comparing two doubles while deciding if the cell is within a sensor
 
HGCalTBCellVertices::HGCalTBCellVertices(){
/*
   ix_ = ix;  
   iv_ = iv;
   sensorsize_ = sensorsize;
*/
   ValidFlag_   = Top.ix_iv_valid(ix_, iv_, sensorsize_);
// Initialize the co-ordinates of a hexagonal cell of side a centred at 0,0 with the top vertex numbered as the first with clockwise increments.
    double x1[] = {0., x_a*a, x_a*a,0., -x_a*a, -x_a*a};
    double y1[] = {a,  y_a*a, -y_a*a, -a, -y_a*a, y_a*a};
    for(int iii=0;iii<6;iii++){// May have to be generalized to deal with polygons of any size
        x_co_FullHex.push_back(x1[iii]);
        y_co_FullHex.push_back(y1[iii]);
       }

  }


std::vector<std::pair<double,double>> HGCalTBCellVertices::GetCellCoordinates(int ix, int iv, int sensorsize){
     ix_ = ix;
     iv_ = iv;
     sensorsize_ = sensorsize;
     ValidFlag_   = Top.ix_iv_valid(ix_, iv_, sensorsize_);
     double vertex_x_tmp=0.,vertex_y_tmp=0.;
     Cell_co.clear();
    if(ValidFlag_){
        for(int iii=0;iii<6;iii++){// May have to be generalized to deal with polygons of any size
           vertex_x_tmp = x_co_FullHex[iii] + ix*x0 + iv*vx0;
           vertex_y_tmp = y_co_FullHex[iii] + iv*vy0;
           if(fabs(vertex_x_tmp) <= Xmax(fabs(vertex_y_tmp)) + delta) Cell_co.push_back(std::make_pair(vertex_x_tmp,vertex_y_tmp));
          }
        return Cell_co; 
       }
     else{
           Cell_co.push_back(std::make_pair(-1.,-1.));
           return Cell_co;
          } 

    }


std::pair<double,double> HGCalTBCellVertices::GetCellCentreCoordinates(int ix, int iv, int sensorsize){
      double centre_x_tmp=0.,centre_y_tmp=0.;  
      ValidFlag_   = Top.ix_iv_valid(ix_, iv_, sensorsize_);
      if(ValidFlag_){
         centre_x_tmp = ix*x0 + iv*vx0;
         centre_y_tmp = iv*vy0;
         return std::make_pair(centre_x_tmp,centre_y_tmp);
        } 
      else return std::make_pair(-1,-1);

     }

double HGCalTBCellVertices::Xmax(double y){
       if(fabs(iv_) <= 3) return 11*x_a*a;
       else return (11*a - y)/(1/(2*x_a));

 } 

/*
void HGCalTBCellVertices::CellType(int ix, int iv, bool ValidFlag){
    bool HalfHex = false;
    x_max=0.;
    y_max=0.;
   }
*/

/*
std::ostream& operator<<(std::ostream& s, const HGCalTBCellVertices& vertices, int type) {

  if(type == 1){
    return s << "This cell is a full hexagon"<<endl;
    for(int iii=0;iii<6;iii++) 
       return s<<" Vertex 1 x co-ordinate = "<<vertices.HX[iii]<<" y co-ordinate = "<<vertices.Hy[iii]<<endl;
    }
  else if(type == 2){
      return s << "This cell is a half hexagon"<<endl;
      for(int iii=0;iii<4;iii++) 
       return s<<" Vertex 1 x co-ordinate = "<<vertices.HX[iii]<<" y co-ordinate = "<<vertices.Hy[iii]<<endl; 
     }  

}
*/
