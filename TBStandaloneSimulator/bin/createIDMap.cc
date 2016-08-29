// --------------------------------------------------------------------------
// File: createIDMap
//
// Description: consolidate various mappings into a single text file
// containing
// sim cell id (cellid), (u, v), (x, y), position id (posid),
// SKIROC number (skiroc), channel number (channel). The text
// HGCal/TBStandaloneSimulator/data/sensor_cellid_uv_map.txt
// is read by HGCCellMap, which provides between various numbers.
//
// Created April 2016 Harrison B. Prosper
// --------------------------------------------------------------------------
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2Poly.h"
#include "TIterator.h"
#include "TList.h"
#include "TText.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/SKIROCParameters.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGeometryConversion.hh"
// --------------------------------------------------------------------------
using namespace std;

namespace
{
int MODEL = 5;               // TB2016 model
int NCELL = 11;              // number of pixels from side to side in sensor
// side length of one pixel (cell) in mm
double CELL_SIDE = 10 * HGCAL_TB_CELL::FULL_CELL_SIDE;
// side length of sensor
double SENSOR_SIDE = NCELL * CELL_SIDE;
double WIDTH = 2 * SENSOR_SIDE; // width of sensor corner to corner
struct Cell {
	int skiroc;
	int channel;
	int sensor_u;
	int sensor_v;
	int celltype;
};
};
// --------------------------------------------------------------------------
int main(int argc, char **argv)
{
	// electronics id to (u, v) mapping file
	string electronicsmapfile("map_FNAL_1234.txt");
	if ( argc > 1 )
		electronicsmapfile = string(argv[1]);

	cout << endl
	     <<  "=> Using (u, v) to (skiroc, channel) map file "
	     << endl
	     <<  "   HGCal/CondObjects/data/"
	     << electronicsmapfile
	     << endl << endl;

	// ----------------------------------------------------------
	// load mapping between (layer, u, v) tp (SKIROC, CHANNEL)
	// ----------------------------------------------------------
	char mfile[1024];
	sprintf(mfile, "$CMSSW_BASE/src/HGCal/CondObjects/data/%s",
	        electronicsmapfile.c_str());
	sprintf(mfile, "%s", gSystem->ExpandPathName(mfile));
	ifstream mapin(mfile);
	if ( ! mapin.good() ) {
		cout <<  "*** cannot open file " << mfile << endl;
		exit(0);
	}
	// dismiss header
	string line;
	getline(mapin, line);
	int skiroc, channel, layer, sensor_u, sensor_v, u, v, celltype;
	map<int, map<pair<int, int>, Cell> > layer2uv2eid;
	vector<int> layers;
	int lastlayer = -1;
	while (mapin
	        >> skiroc >> channel >> layer
	        >> sensor_u >> sensor_v
	        >> u >> v
	        >> celltype) {
		Cell cell;
		cell.skiroc   = skiroc;
		cell.channel  = channel;
		cell.sensor_u = sensor_u;
		cell.sensor_v = sensor_v;
		cell.celltype = celltype;
		pair<int, int> uv(u, v);

		if ( layer2uv2eid.find(layer) == layer2uv2eid.end() )
			layer2uv2eid[layer] = map<pair<int, int>, Cell>();
		layer2uv2eid[layer][uv] = cell;

		if ( layer != lastlayer ) layers.push_back(layer);
		lastlayer = layer;
	}

	// ----------------------------------------------------------
	// create a 2-D histogram with hexagonal bins, a
	// subset of which lie within the hexagonal boundary
	// that defines a sensor
	// ----------------------------------------------------------
	HGCSSGeometryConversion geom(MODEL, WIDTH, CELL_SIDE);
	geom.initialiseHoneyComb(WIDTH, CELL_SIDE);
	TH2Poly* map = geom.hexagonMap();

	// ----------------------------------------------------------
	// create a single hexagonal bin to represent sensor.
	// we shall use this hexagon region to determine which
	// sim cellids lie within a sensor.
	// ----------------------------------------------------------
	TH2Poly hsensor;
	hsensor.SetName("hsensor");
	hsensor.SetTitle("sensor");

	// make slightly smaller so that we identify
	// mouse bitten cells
	double S = 0.98 * SENSOR_SIDE;
	double H = S * sqrt(3) / 2; // center to side distance
	double x[7], y[7];
	x[0] = -S / 2;
	y[0] = -H;
	x[1] = -S;
	y[1] =  0;
	x[2] = -S / 2;
	y[2] =  H;
	x[3] =  S / 2;
	y[3] =  H;
	x[4] =  S;
	y[4] =  0;
	x[5] =  S / 2;
	y[5] = -H;
	x[6] = -S / 2;
	y[6] = -H;
	hsensor.AddBin(7, x, y);

	//get the single hexagonal bin that represents
	//the boundary of sensor
	TIter it(hsensor.GetBins());
	TH2PolyBin* sensor = (TH2PolyBin*)it();

	// loop over bins in 2-D histogram, determine which
	// ones lie within sensor, and write out the bin
	// numbers (cellid) and the (x,y) centers of each
	// pixel.
	gStyle->SetPalette(1);
	gStyle->SetOptStat("");

	TCanvas csensor("sensor_cellids", "cellid", 600, 600);
	map->GetXaxis()->CenterTitle();
	map->GetXaxis()->SetTitle("#font[12]{x} axis");
	map->GetYaxis()->CenterTitle();
	map->GetYaxis()->SetTitle("#font[12]{y} axis");

	// make a sensor with the correct size (for plotting)
	TH2Poly hsensortrue;
	hsensortrue.SetName("hsensortrue");
	S = SENSOR_SIDE;
	H = S * sqrt(3) / 2; // center to side distance
	x[0] = -S / 2;
	y[0] = -H;
	x[1] = -S;
	y[1] =  0;
	x[2] = -S / 2;
	y[2] =  H;
	x[3] =  S / 2;
	y[3] =  H;
	x[4] =  S;
	y[4] =  0;
	x[5] =  S / 2;
	y[5] = -H;
	x[6] = -S / 2;
	y[6] = -H;
	hsensortrue.AddBin(7, x, y);
	hsensortrue.SetMinimum(0.0);
	hsensortrue.SetMaximum(1.0);
	hsensortrue.SetBinContent(1, 0.7);

	csensor.cd();
	map->Draw();
	hsensortrue.Draw("col same");
	map->Draw("same");

	// (u, v) plot
	TCanvas cuv("sensor_u_v", "u, v", 600, 600);
	hsensortrue.SetBinContent(1, 0.0);
	hsensortrue.Draw();
	hsensortrue.SetBinContent(1, 0.7);

	cuv.cd();
	map->Draw();
	hsensortrue.Draw("col same");
	map->Draw("same");

	// ----------------------------------------------------------
	// loop over cells
	// ----------------------------------------------------------
	char record[80];
	ofstream sout("sensor_map.txt");
	sprintf(record,
	        "%5s %5s %5s %5s %5s %5s %5s "
	        "%5s %5s %5s "
	        "%9s %9s",
	        "layer", "posid", "cid", "s_u", "s_v", "u", "v",
	        "ski", "chan", "ctype",
	        "x", "y");
	cout << record << endl;
	sout << record << endl;

	TList* bins = map->GetBins();
	TText text;
	text.SetTextSize(0.02);
	text.SetTextAlign(22);  // centered

	HGCalTBTopology topology;
	HGCalTBCellVertices vertices;

	// in offline code, layers start at 1
	// hard code layer and channel count for now
	int ncells = 128;
	for(size_t c = 0; c < layers.size(); c++) {
		TIter next(bins);
		int layer = layers[c];

		// number of cells in y (either 12 or 11)
		bool new_column = true;
		bool odd_column = true;
		int colnumber = 0;
		int u_start = 1;
		int v_start = 8;
		int number = 1;
		u = 0;
		v = 0;
		while ( TH2PolyBin* bin = (TH2PolyBin*)next() ) {
			// We get the starting
			// values of (u, v) per column as follows:
			//   1. every two columns, increment u
			//   2. every column, decrement v
			// Thereafter:
			//   1. decrement u

			int binnumber = bin->GetBinNumber();
			new_column = binnumber == number;
			if ( new_column ) {
				// decrement v
				v_start--;

				colnumber++;
				if ( colnumber % 2 == 1 ) u_start++;

				// initialize (u, v) for current column
				u = u_start + 1;
				v = v_start;

				new_column = false;
				if ( odd_column )
					number += 12;
				else
					number += 11;
				odd_column = !odd_column;
			}

			// decrement u
			u--;

			// map to skiroc number and channel
			pair<int, int> uv(u, v);
			Cell cell = layer2uv2eid[layer][uv];

			if ( ! topology.iu_iv_valid(layer,
			                            cell.sensor_u,
			                            cell.sensor_v,
			                            u, v,
			                            ncells) ) continue;

			pair<double, double> pos
			    = vertices.GetCellCentreCoordinates(layer,
			                                        cell.sensor_u,
			                                        cell.sensor_v,
			                                        u, v,
			                                        ncells);
			double x = pos.first;
			double y = pos.second;
			x *= 10; // change to mm
			y *= 10;

			// posid is the position identifier, which specifies whether the
			// cell is an interior cell or a boundary cell.
			int posid = 0;
			if ( ! sensor->IsInside(x, y) ) {
				// boundary cell, determine type
				if ( v > 3 ) {
					// either lower or upper left
					if ( u < -3 )
						posid = 3; // upper
					else
						posid = 4; // lower
				} else if ( v < -3 ) {
					// either lower or upper right
					if ( u > 3 )
						posid = 6; // lower
					else
						posid = 1; // upper
				} else {
					// either lower or upper
					if ( u > 0 )
						posid = 5; // lower
					else
						posid = 2; // upper
				}
			}

			csensor.cd();
			sprintf(record, "%d", binnumber);
			text.DrawText(x, y, record);
			cuv.cd();
			sprintf(record, "%d,%d", u, v);
			text.DrawText(x, y, record);

			sprintf(record,
			        "%5d %5d %5d %5d %5d %5d %5d "
			        "%5d %5d %5d "
			        "%9.3f %9.3f",
			        layer, posid, binnumber,
			        cell.sensor_u, cell.sensor_v,
			        u, v,
			        cell.skiroc, cell.channel, cell.celltype,
			        x, y);
			cout << record << endl;
			sout << record << endl;
		}
	}
	sout.close();

	map->SetTitle("TB2016 Standalone Simulator Cell IDs");
	csensor.Update();
	csensor.SaveAs(".png");

	map->SetTitle("TB2016 Sensor (u,v) Coordinates");
	cuv.Update();
	cuv.SaveAs(".png");

	return 0;
}
