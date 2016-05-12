#include "HGCal/TBStandaloneSimulator/interface/HGCSSGenParticle.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>


void HGCSSGenParticle::Print(std::ostream & aOs) const
{
	aOs << std::setprecision(6)
	    << "====================================" << std::endl
	    << " = time " << time_ << " ns" << std::endl
	    << " = position " << xpos_ << " " << ypos_ << " " << zpos_ << " mm" << std::endl
	    << " = Mass " << mass_ << " MeV" << std::endl
	    << " = momentum " << px_ << " " <<  py_ << " " << pz_ << " MeV" << std::endl
	    << " = pdgid " << pdgid_ << std::endl
	    << " = charge " << charge_ << std::endl
	    << " = G4trackID " << trackID_ << std::endl
	    << "====================================" << std::endl;

}
void HGCSSGenParticle::Print(const unsigned idx,
                             std::ostream & aOs) const
{
	aOs << std::setprecision(6)
	    << "========= GenParticle " << idx << " =========" << std::endl;
	Print(aOs);

}

