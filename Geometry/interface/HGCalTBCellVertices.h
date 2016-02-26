#ifndef HGCAL_GEOMETRY_HGCALTBCELLVERTICES_H
#define HGCAL_GEOMETRY_HGCALTBCELLVERTICES_H

/** \class HGCalTBCellVertices HGCal/Geometry/interface/HGCalTBCellVertices.h HGCal/Geometry/src/HGCalTBCellVertices.cc

\brief definition of a cell

This class implements the elementary cell in HGCal
*/

#include <cmath>
#include <vector>

class HGCalTBCellVertices
{
public:

	HGCalTBCellVertices(); ///<  Constructor from cell ix & iv, valid sensorSizes are 128 and 256
	std::vector<std::pair<double, double>> GetCellCoordinates(int ix, int iv, int sensorsize) const;
	std::pair<double, double> GetCellCentreCoordinates(int ix, int iv, int sensorsize);
//  void CellType(int ix, int iv, bool ValidFlag);// 1 for full hex, 2 for half hex and 3 for the pentagons(to be implemented later)
private:
	const double a = 1; ///< Size in terms of 1 unit of x/y co-ordinates of a cell side
	const double x_a = sqrt(3) / 2; // cosine pi/6
	const double y_a = 1 / 2.; // sine pi/6
	const double vy_a = 3. / 2;

	std::vector<double> x_co_FullHex, y_co_FullHex; // stores the initial x,y coordinates of a hexagonal cell
	std::vector<std::pair<double, double>> Cell_co;
// Translation in x,v co-ordinates in terms of cartesian x,y.
	double  x0 = 2 * x_a * a; //Translation in Cartesian x for 1 unit of ix
	double vx0 = x_a * a; // Cartesian x component of translation for 1 unit of iv
	double vy0 = vy_a * a; // Cartesian y component of translation for 1 unit of iv
	int ix_, iv_, sensorsize_;
	bool ValidFlag_;
	double x_max, y_max;

	double Xmax(double y);// returns the max x value for a cell to be in the given sensor

};

//std::ostream& operator<<(std::ostream&,const HGCalTBCellVertices& Cell_Ix_Iv);

#endif

