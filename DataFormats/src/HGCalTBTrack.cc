#include "HGCal/DataFormats/interface/HGCalTBTrack.h"
#include <cassert>

HGCalTBTrack::HGCalTBTrack(void)
{
}

HGCalTBTrack::HGCalTBTrack(const float *v){
//	assert(size == getSize()); ///\todo to be turned into exception
	
	chi2	= v[0];
	x0		= v[1];
	y0		= v[2];
	m_x		= v[3];
	m_y		= v[4];
	m_x_err = v[5];
	m_y_err = v[6];

	return;

}

std::vector<float> HGCalTBTrack::getRaw(void)
{
	std::vector<float> v;

	v.push_back(chi2);
	v.push_back(x0);
	v.push_back(y0);
	v.push_back(m_x);
	v.push_back(m_y);
	v.push_back(m_x_err);
	v.push_back(m_y_err);

	return v;
}
