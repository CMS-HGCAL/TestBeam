#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "math.h"
#include <stdlib.h>
#include <iostream>
#define PI 3.14159265

bool HGCalTBTopology::iu_iv_valid(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorSize) const
{
	bool Is_Valid_sensor_iu_iv = true;
	if( abs(sensor_iu+sensor_iv)>1 || abs(sensor_iu)>1 || abs(sensor_iv)>1 )
	  Is_Valid_sensor_iu_iv=false; //!!ASSUMES NO MORE THAN 7 SEVEN MODULES PER LAYER!!
        int aiv = abs(iv);
        int iuc = (iv < 0) ? (-iu) : (iu);
	if(Is_Valid_sensor_iu_iv) {
		if(sensorSize == 128) {
			if (iv == 0) return (iu >= -5 && iu <= 5);
			else if (aiv == 1) return (iuc >= -6 && iuc <= 5);
			else if (aiv == 2) return (iuc >= -6 && iuc <= 4);
			else if (aiv == 3) return (iuc >= -7 && iuc <= 4);
			else if (aiv == 4) return (iuc >= -7 && iuc <= 3);
			else if (aiv == 5) return (iuc >= -6 && iuc <= 1);
			else if (aiv == 6) return (iuc >= -5 && iuc <= -1);
			else if (aiv == 7) return (iuc == -3 || iuc == -4);
			else return false;
		} else return false;
	} else return false;
}

double HGCalTBTopology::Cell_Area(int cell_type) const
{
	double a = HGCAL_TB_CELL::FULL_CELL_SIDE;
	double b = HGCAL_TB_CELL::CALIB_PAD_SIDE;
	double mouse_bite = HGCAL_TB_CELL::MOUSE_BITE_SIDE;

	double area_full_hex = 0.5 * 3 * sqrt(3) * pow(a, 2);
	double area_calib_pad = 0.5 * 3 * sqrt(3) * pow(b, 2);
	double area_half_hex = area_full_hex / 2;
	double area_mouse_bite = area_half_hex - (0.5 * pow((a - mouse_bite), 2) / tan(PI / 6));
	double area_outer_calib_pad = area_full_hex - area_calib_pad;
	double area_merged_cell = area_full_hex + area_half_hex;

	if (cell_type == 0) return area_full_hex;
	else if (cell_type == 1) return area_calib_pad;
	else if (cell_type == 2) return area_half_hex;
	else if (cell_type == 3) return area_mouse_bite;
	else if (cell_type == 4) return area_outer_calib_pad;
	else if (cell_type == 5) return area_merged_cell;
	else return -1.; //signifies an invalid cell type
}


std::set<HGCalTBDetId> HGCalTBTopology::getNeighboringCellsDetID(HGCalTBDetId detid, int sensorSize, int maxDistance, const HGCalElectronicsMap& emap) const
{
	int layer=detid.layer();
	std::set<HGCalTBDetId> detids;
	for(int u=-maxDistance; u<=maxDistance; u++){
		for(int v=-maxDistance; v<=maxDistance; v++){
			if( (u==0 && v==0) || abs(u+v)>maxDistance ) continue;
			int iU=detid.sensorIU();
			int iV=detid.sensorIV();
			int iu=detid.iu()+u;
			int iv=detid.iv()+v;
			if( iu_iv_valid( layer, iU, iV, iu, iv, sensorSize)!=true )
				continue;
			for(unsigned int cellType=0; cellType<6; cellType++){
				HGCalTBDetId did( layer, iU, iV, iu, iv, cellType );
				if( emap.existsDetId(did)==true ){
					detids.insert(did);
					if( cellType==1 ){
						HGCalTBDetId didbis( layer, iU, iV, iu, iv, 4);
						detids.insert(didbis);
					}
					break;
				}
			}
		}
	}
	return detids;
}

