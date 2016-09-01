#ifndef HGCCELLMAP_H
#define HGCCELLMAP_H
/** \class

	\author Harrison B. Prosper
	map from sim cell id to (u,v) coordinates
	and from (u,v) to (x, y) coordinates assuming
	cell is centered at (0, 0).
 */
#include <map>
#include <string>
#include <vector>
#include "TH2Poly.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGeometryConversion.hh"

class HGCCellMap
{
public:
	struct Cell {
		int skiroc;
		int channel;
		int celltype;
		int posid;
		int u;
		int v;
		double x;
		double y;
		double z;
		double count;
		bool operator<(HGCCellMap::Cell& o)
		{
			// SKIROC 2 occurs before SKIROC 1
			int ski   = skiroc % 2;
			int o_ski = o.skiroc % 2;
			int lhs   = 64 * ski + channel;
			int rhs   = 64 * o_ski + o.channel;
			return lhs < rhs;
		}
	};
	HGCCellMap(std::string inputFilename = "");
	virtual ~HGCCellMap();

	std::pair<int, int>       operator()(size_t cellid);
	std::pair<double, double> operator()(std::pair<int, int>& uv);
	std::pair<double, double> uv2xy(int u, int v);
	std::pair<int, int>       xy2uv(double x, double y);
	int                       posid(int u, int v);

	std::pair<int, int>       uv2eid(int layer,
	                                 int sensor_u, int sensor_v,
	                                 int u, int v);
	int                       celltype(int layer,
	                                   int sensor_u, int sensor_v,
	                                   int u, int v);

	std::vector<HGCCellMap::Cell> cells(int layer = 1,
	                                    int sensor_u = 0,
	                                    int sensor_v = 0);
	TH2Poly*                  honeycomb()
	{
		return _map;
	}
private:
	std::map<size_t, std::pair<int, int> > _uvmap;
	std::map<std::pair<int, int>, int> _posid;
	std::map<std::pair<int, int>, std::pair<double, double> > _xymap;

	std::map<uint32_t, int> _celltype;
	std::map<uint32_t, std::pair<int, int> > _eidmap;
	std::map<uint32_t, std::vector<HGCCellMap::Cell> > _cells;

	HGCSSGeometryConversion _geom;
	TH2Poly* _map;
};

#endif
