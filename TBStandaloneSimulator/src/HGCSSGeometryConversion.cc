#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "TH2Poly.h"
#include "TMath.h"

#include "HGCal/TBStandaloneSimulator/interface/HGCSSGeometryConversion.hh"


HGCSSGeometryConversion::HGCSSGeometryConversion(const int model,
        const double width,
        const double cellsize)
	: model_(model),
	  width_(width),
	  cellsize_(cellsize)
{
}

HGCSSGeometryConversion::~HGCSSGeometryConversion()
{
	hexaGeom.clear();
	squareGeom.clear();
}

void HGCSSGeometryConversion::initialiseSquareMap(const double xymin,
        const double side)
{
	initialiseSquareMap(squareMap(), xymin, side, true);
	fillXY(squareMap(), squareGeom);
}

void HGCSSGeometryConversion::initialiseSquareMap(TH2Poly *map,
        const double xymin,
        const double side,
        bool print)
{
	unsigned nx = static_cast<unsigned>(xymin * 2. / side);
	unsigned ny = nx;
	unsigned i, j;
	Double_t x1, y1, x2, y2;
	Double_t dx = side, dy = side;
	x1 = -1.*xymin;
	x2 = x1 + dx;

	for (i = 0; i < nx; i++) {
		y1 = -1.*xymin;
		y2 = y1 + dy;
		for (j = 0; j < ny; j++) {
			map->AddBin(x1, y1, x2, y2);
			y1 = y2;
			y2 = y1 + dy;
		}
		x1 = x2;
		x2 = x1 + dx;
	}

	if (print) {
		std::cout <<  " -- Initialising squareMap with parameters: " << std::endl
		          << " ---- xymin = " << -1.*xymin << ", side = " << side
		          << ", nx = " << nx << ", ny=" << ny
		          << std::endl;
	}

}

void HGCSSGeometryConversion::initialiseHoneyComb(const double width,
        const double side)
{
	initialiseHoneyComb(hexagonMap(), width, side, true);
	fillXY(hexagonMap(), hexaGeom);
}

void HGCSSGeometryConversion::initialiseHoneyComb(TH2Poly *map,
        const double width,
        const double side,
        bool print)
{
	// Center a cell at (x,y)=(0,0) and ensure
	// coverage up to/past width/2 in all 4 directions,
	// assuming each cell is lying on a side.

	unsigned ncellwide = 11;
	unsigned ny = ncellwide + 1;
	unsigned nx = ncellwide + 4;
	double xstart = -((double)ncellwide + 0.5) * side;
	double ystart = -((double)ncellwide + 1) * side * sqrt(3) / 2;
	if (print) {
		std::cout << " -- Initialising HoneyComb with parameters: " << std::endl
		          << " ---- (xstart,ystart) = (" << xstart << "," << ystart << ")"
		          << ", side = " << side << ", nx = " << nx << ", ny=" << ny << std::endl;
	}
	myHoneycomb(map, xstart, ystart, side, ny, nx);
}

////////////////////////////////////////////////////////////////////////////////
/// Copied and modified from TH2Poly.cxx
/// Bins the histogram using a honeycomb structure
/// 90 degree rotation, side up instead of corner up

void HGCSSGeometryConversion::myHoneycomb(TH2Poly* map,
        Double_t xstart, Double_t ystart,
        Double_t a,  // side length
        Int_t k,     // # hexagons in a column
        Int_t s)     // # columns
{
	// Add the bins
	Double_t numberOfHexagonsInAColumn;
	Double_t x[6], y[6];
	Double_t xloop, yloop, ytemp;
	xloop = xstart;
	yloop = ystart + a * TMath::Sqrt(3) / 2.0;
	for (int sCounter = 0; sCounter < s; sCounter++) {

		ytemp = yloop; // Resets the temp variable

		// Determine the number of hexagons in that column
		if(sCounter % 2 == 0) {
			numberOfHexagonsInAColumn = k;
		} else {
			numberOfHexagonsInAColumn = k - 1;
		}

		for (int kCounter = 0; kCounter <  numberOfHexagonsInAColumn; kCounter++) {

			// Go around the hexagon
			x[0] = xloop;
			y[0] = ytemp;
			x[1] = x[0] + a / 2.0;
			y[1] = y[0] + a * TMath::Sqrt(3) / 2.0;
			x[2] = x[1] + a;
			y[2] = y[1];
			x[3] = x[2] + a / 2.0;
			y[3] = y[1] - a * TMath::Sqrt(3) / 2.0;;
			x[4] = x[2];
			y[4] = y[3] - a * TMath::Sqrt(3) / 2.0;;
			x[5] = x[1];
			y[5] = y[4];

			map->AddBin(6, x, y);

			// Go up
			ytemp += a * TMath::Sqrt(3);
		}

		// Increment the starting position
		if (sCounter % 2 == 0) yloop += a * TMath::Sqrt(3) / 2.0;
		else                 yloop -= a * TMath::Sqrt(3) / 2.0;
		xloop += 1.5 * a;
	}
}


void HGCSSGeometryConversion::fillXY(TH2Poly* hist,
                                     std::map<int, std::pair<double, double> >&
                                     geom)
{
	TIter next(hist->GetBins());
	TObject *obj = 0;
	TH2PolyBin *polyBin = 0;
	geom.clear();

	while ((obj = next())) {
		polyBin = (TH2PolyBin*)obj;
		int id = polyBin->GetBinNumber();
		std::pair<double, double> xy = std::pair<double, double>((polyBin->GetXMax() + polyBin->GetXMin()) / 2., (polyBin->GetYMax() + polyBin->GetYMin()) / 2.);
		geom.insert(std::pair<unsigned, std::pair<double, double> >(id, xy));
	}
	std::cout << " -- Check geomMap: size = " << geom.size() << std::endl;
}






