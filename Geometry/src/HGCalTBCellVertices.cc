#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include <stdlib.h>
#include "math.h"
#include <iostream>
#include <iomanip>
#define PI 3.14159265
#define TEST_BEAM_LAYER_ROTATION -PI/2
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


std::vector<std::pair<double, double>> HGCalTBCellVertices::GetCellCoordinates(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorsize, bool flipX)
{
	bool ValidFlag   = Top.iu_iv_valid(layer, sensor_iu, sensor_iv, iu, iv, sensorsize);
	double vertex_x_tmp = 0., vertex_y_tmp = 0.;
	Cell_co.clear();
	if(ValidFlag) {
		for(int iii = 0; iii < 6; iii++) { // May have to be generalized to deal with polygons of any size
			vertex_x_tmp = x_co_FullHex[iii] + iu * x0 + iv * vx0;
			vertex_y_tmp = y_co_FullHex[iii] + iv * vy0;
//The general strategy is to translate starting from the central hexagonal cell to the iu,iv desired. If any vertex goes out of the sensor boundary its cordinates are not filled into the vector of pairs.
			if(fabs(vertex_x_tmp) <= Xmax(iv, fabs(vertex_y_tmp)) + delta) {
				vertex_x_tmp += sensor_iu*X0 + sensor_iv*VX0;
				vertex_y_tmp += sensor_iv*VY0;
				auto point = RotateLayer(std::make_pair(vertex_x_tmp, vertex_y_tmp), TEST_BEAM_LAYER_ROTATION);
//				if(flipX==true) point.first=-point.first;
				Cell_co.push_back(point);
			}
		}
		return Cell_co;
	} else {
		Cell_co.push_back(std::make_pair(-123456, -123456)); //iu_iv_Valid() is sufficient to decide if a given iu,iv is a valid sensor index but this is done if some future need may arise.
		return Cell_co;
	}

}


std::pair<double, double> HGCalTBCellVertices::GetCellCentreCoordinates(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorsize, bool flipX)
{
	double centre_x_tmp = 0., centre_y_tmp = 0.;
	bool ValidFlag   = Top.iu_iv_valid(layer, sensor_iu, sensor_iv, iu, iv, sensorsize);
	if(ValidFlag) {    
		centre_x_tmp = ( (iu * x0 + iv * vx0) < 0) ? ((iu * x0 + iv * vx0) + delta) : ((iu * x0 + iv * vx0) - delta)  ;
		centre_y_tmp = ( (iv * vy0) < 0) ? ((iv * vy0) + delta) : ((iv * vy0) - delta);
		
		centre_x_tmp += sensor_iu*X0 + sensor_iv*VX0;
		centre_y_tmp += sensor_iv*VY0;
		auto point = RotateLayer(std::make_pair(centre_x_tmp, centre_y_tmp), TEST_BEAM_LAYER_ROTATION);
//		if(flipX==true) point.first = - point.first;
		return point;
    
	} else return std::make_pair(-123456, -123456); //iu_iv_Valid() is sufficient to decide if a given iu,iv is a valid sensor index but this is done if some future need may arise.

}

std::pair<int, int> HGCalTBCellVertices::GetCellIUIVCoordinates(double x, double y) { 
  double iv = x * 1/vy0;
  double iu = (y-vx0*iv)/x0;

  return std::make_pair(round(iu), round(iv));
  
}


std::pair<int, int> HGCalTBCellVertices::GetSensorIUIVCoordinates(double X, double Y) {  
  double iV = X * 1/VY0;
  double iU = (Y-VX0*iV)/X0;

  return std::make_pair(round(iU), round(iV));
  
}


double HGCalTBCellVertices::Xmax(int iv, double y)
{
	if(fabs(iv) <= 3) return 11 * x_a * a;
	else return (11 * a - y) / (1 / (2 * x_a));
}

std::pair<double, double> HGCalTBCellVertices::RotateLayer(std::pair<double, double> Vertex, double Angle)
{
	double X_new = (Vertex.first) * cos(Angle) - (Vertex.second) * sin(Angle);
	double Y_new = (Vertex.first) * sin(Angle) + (Vertex.second) * cos(Angle);
	return std::make_pair(-X_new, Y_new);// The negative sign for the x coordinate is to account for the TB cartesian coordinate system.
}

// To be added if for reconstruction it is useful to simply know if a cell is full hex, half hex or mouse-bitten
/*
void HGCalTBCellVertices::CellType(int iu, int iv, bool ValidFlag){
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
