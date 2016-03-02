#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include <stdlib.h>

bool HGCalTBTopology::ix_iv_valid(int layer, int sensor_ix, int sensor_iv, int ix, int iv, int sensorSize) const
{
	int aiv = abs(iv);
	int ixc = (iv < 0) ? (-ix) : (ix);
             if(layer <= 28 && sensor_ix == 0 && sensor_iv == 0){
                 if(sensorSize == 128) {
	            if (iv == 0) return (ix >= -5 && ix <= 5);
		    else if (aiv == 1) return (ixc >= -6 && ixc <= 5);
		    else if (aiv == 2) return (ixc >= -6 && ixc <= 4);
	            else if (aiv == 3) return (ixc >= -7 && ixc <= 4);
	            else if (aiv == 4) return (ixc >= -7 && ixc <= 3);
		    else if (aiv == 5) return (ixc >= -6 && ixc <= 1);
		    else if (aiv == 6) return (ixc >= -5 && ixc <= -1);
		    else if (aiv == 7) return (ixc == -3 || ixc == -4);
		    else return false;
                   }
                  else return false;
                }
              else return false;
}
