#ifndef _g4sihit_hh_
#define _g4sihit_hh_

#include <vector>
class G4SiHit
{

public:
	G4SiHit();
	~G4SiHit();

	double energy;
	double time;
	unsigned layer;
	int pdgId;
	double hit_x;
	double hit_y;
	double hit_z;
	int trackId;
	int parentId;

private:

};


typedef std::vector<G4SiHit> G4SiHitVec;



#endif
