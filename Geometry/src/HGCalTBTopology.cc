#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include <stdlib.h>

bool HGCalTBTopology::iu_iv_valid(int layer, int sensor_iu, int sensor_iv, int iu, int iv, int sensorSize) const
{
	int aiv = abs(iv);
	int iuc = (iv < 0) ? (-iu) : (iu);
	if(layer <= 28 && sensor_iu == 0 && sensor_iv == 0) {
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
