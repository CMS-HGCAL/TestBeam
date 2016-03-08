#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include <stdlib.h>
#include "math.h"
#include <iostream>
#include <iomanip>
#define PI 3.14159265
using namespace std;

HGCalTBTopology Top;
double delta = 0.00001;//needed for comparing two doubles while deciding if the cell is within a sensor

HGCalTBCellVertices::HGCalTBCellVertices()
{
// Initialize the co-ordinates of a hexagonal cell of side a centred at 0,0 with the top vertex numbered as the first with clockwise increments.
	double x1[] = {0., x_a * a, x_a * a, 0., -x_a * a, -x_a * a};
	double y1[] = {a,  y_a * a, -y_a * a, -a, -y_a * a, y_a * a};
	for(int iii = 0; iii < 6; iii++) { // May have to be generalized to deal with polygons of any size
		x_co_FullHex.push_back(x1[iii]);
		y_co_FullHex.push_back(y1[iii]);
	}

}


std::vector<std::pair<double, double>> HGCalTBCellVertices::GetCellCoordinates(int layer, int sensor_ix, int sensor_iv, int ix, int iv, int sensorsize)
{
	bool ValidFlag   = Top.ix_iv_valid(layer, sensor_ix, sensor_iv, ix, iv, sensorsize);
	double vertex_x_tmp = 0., vertex_y_tmp = 0.;
	Cell_co.clear();
	if(ValidFlag) {
		for(int iii = 0; iii < 6; iii++) { // May have to be generalized to deal with polygons of any size
			vertex_x_tmp = x_co_FullHex[iii] + ix * x0 + iv * vx0;
			vertex_y_tmp = y_co_FullHex[iii] + iv * vy0;
//The general strategy is to translate starting from the central hexagonal cell to the ix,iv desired. If any vertex goes out of the sensor boundary its cordinates are not filled into the vector of pairs.
			if(fabs(vertex_x_tmp) <= Xmax(iv, fabs(vertex_y_tmp)) + delta) Cell_co.push_back(std::make_pair(vertex_x_tmp, vertex_y_tmp));
		}
		return Cell_co;
	} else {
		Cell_co.push_back(std::make_pair(-123456, -123456)); //ix_iv_Valid() is sufficient to decide if a given ix,iv is a valid sensor index but this is done if some future need may arise.
		return Cell_co;
	}

}


std::pair<double, double> HGCalTBCellVertices::GetCellCentreCoordinates(int layer, int sensor_ix, int sensor_iv, int ix, int iv, int sensorsize)
{
	double centre_x_tmp = 0., centre_y_tmp = 0.;
	bool ValidFlag   = Top.ix_iv_valid(layer, sensor_ix, sensor_iv, ix, iv, sensorsize);
	if(ValidFlag) {
		centre_x_tmp = ix * x0 + iv * vx0;
		centre_y_tmp = iv * vy0;
		return std::make_pair(centre_x_tmp, centre_y_tmp);
	} else return std::make_pair(-123456, -123456); //ix_iv_Valid() is sufficient to decide if a given ix,iv is a valid sensor index but this is done if some future need may arise.

}

double HGCalTBCellVertices::Xmax(int iv, double y)
{
	if(fabs(iv) <= 3) return 11 * x_a * a;
	else return (11 * a - y) / (1 / (2 * x_a));
}

// To be added if for reconstruction it is useful to simply know if a cell is full hex, half hex or mouse-bitten
/*
void HGCalTBCellVertices::CellType(int ix, int iv, bool ValidFlag){
    bool HalfHex = false;
    x_max=0.;
    y_max=0.;
   }
*/

// To be added after finalizing what all we need to see printed out.
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
