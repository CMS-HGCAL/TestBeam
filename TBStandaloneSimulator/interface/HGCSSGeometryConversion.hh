#ifndef HGCSSGeometryConversion_h
#define HGCSSGeometryConversion_h

#include <map>
#include "TH2Poly.h"

class HGCSSGeometryConversion
{
public:
	HGCSSGeometryConversion() {}
	HGCSSGeometryConversion(const int model,
	                        const double width,
	                        const double cellsize);

	~HGCSSGeometryConversion();

	TH2Poly *hexagonMap()
	{
		static TH2Poly hc;
		return &hc;
	};

	TH2Poly *squareMap()
	{
		static TH2Poly hsq;
		return &hsq;
	};

	std::map<int, std::pair<double, double> > hexaGeom;
	std::map<int, std::pair<double, double> > squareGeom;

	void initialiseSquareMap(const double xymin, const double side);

	void initialiseSquareMap(TH2Poly* hmap,
	                         const double xymin,
	                         const double side, bool print);

	void initialiseHoneyComb(const double xymin, const double side);

	void initialiseHoneyComb(TH2Poly* hmap,
	                         const double xymin,
	                         const double side, bool print);

	void fillXY(TH2Poly* hist, std::map<int, std::pair<double, double> > & geom);

	inline double getXYwidth() const
	{
		return width_;
	};

	inline void setXYwidth(double width)
	{
		width_ = width;
	};

	inline double cellSize() const
	{
		return cellsize_;
	};

private:

	void myHoneycomb(TH2Poly* pmap,
	                 Double_t xstart,
	                 Double_t ystart,
	                 Double_t a,  // side length
	                 Int_t k,     // # hexagons in a column
	                 Int_t s);    // # columns
	unsigned model_;
	double width_;
	double cellsize_;
};



#endif
