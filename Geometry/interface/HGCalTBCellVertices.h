#ifndef HGCAL_GEOMETRY_HGCALTBCELLVERTICES_H
#define HGCAL_GEOMETRY_HGCALTBCELLVERTICES_H

#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "math.h"
#include <vector>

/**
 * \class HGCal/Geometry/interface/HGCalTBCellVertices.h
 *
 * \brief This class implements the local coordinate system
 *
 * The local coordinate system is a 2D non-hortogonal system. The axes are \b iu and \b iv.
 * The orientation of the axis is \b iu
 *
 * \author Rajdeep Mohan Chatterjee, Shervin Nourbakhsh
 *
 */

class HGCalTBCellVertices
{
public:

	HGCalTBCellVertices();	///< Constructor from cell \b iu & \b iv, valid sensorSizes are 128 and 256

	std::vector<std::pair<double, double>> GetCellCoordinates(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorsize, bool flipX = false); ///< returns the coordinates of each vertex of cell in the lab frame \b (x,y)

	inline std::vector<std::pair<double, double>> GetCellCoordinatesForPlots(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorsize)
	{
		return GetCellCoordinates(layer, sensor_iu, sensor_iv, iu, iv, sensorsize, true);
	};  ///< returns the coordinates of each vertex of cell in the lab frame \b (x,y)

	std::pair<double, double> GetCellCentreCoordinates(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorsize, bool flipX = false); ///< returns the center of the cell in absolute coordinates: \b (x,y)
	std::pair<int, int> GetCellIUIVCoordinates(double x, double y);
	std::pair<int, int> GetSensorIUIVCoordinates(double X, double Y);

	inline 	std::pair<double, double> GetCellCentreCoordinatesForPlots(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorsize)
	{
		return GetCellCentreCoordinates(layer, sensor_iu, sensor_iv, iu, iv, sensorsize, true);
	}; ///< returns the center of the cell in absolute coordinates: \b (x,y)


//  void CellType(int iu, int v, bool ValidFlag);// 1 for full hex, 2 for half hex and 3 for the pentagons(to be implemented later)
private:
	double a = HGCAL_TB_CELL::FULL_CELL_SIDE; // Size in terms of 1 unit of x/y co-ordinates of a cell side which is 0.064 cm
	double A = 11*a; // One side of a full sensor(neglecting the cut at the MB)
	double x_a = sqrt(3) / 2; // cosine pi/6
	double y_a = 1 / 2.; // sine pi/6
	double vy_a = 3. / 2;

	std::vector<double> x_co_FullHex, y_co_FullHex; // stores the initial x,y coordinates of a hexagonal cell
	std::vector<std::pair<double, double>> Cell_co;
// Translation in u,v co-ordinates in terms of TB cartesian -x,y.
	double  x0 = 2 * x_a * a; //Translation in Cartesian x for 1 unit of iu
	double vx0 = x_a * a; // Cartesian x component of translation for 1 unit of iv
	double vy0 = vy_a * a; // Cartesian y component of translation for 1 unit of iv
// Translation in Sensor_u, Sensor_v co-ordinates in terms of TB cartesian -x,y.
        double  X0 = 2 * x_a * A; //Translation in Cartesian x for 1 unit of Sensor_iu
        double VX0 = x_a * A; // Cartesian x component of translation for 1 unit of Sensor_iv
        double VY0 = vy_a * A; // Cartesian y component of translation for 1 unit of Sensor_iv

	double Xmax(int iv, double y);// returns the max x value for a cell to be in the given sensor

	std::pair<double, double> RotateLayer(std::pair<double, double>, double Angle);

};

//std::ostream& operator<<(std::ostream&,const HGCalTBCellVertices& Cell_U_V);

#endif
