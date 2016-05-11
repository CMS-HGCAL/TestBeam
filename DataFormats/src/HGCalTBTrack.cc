#include "HGCal/DataFormats/interface/HGCalTBTrack.h"
#include <cassert>

double MicronToCM    = 0.0001;// Tracker information has microns as the length unit, conversion factor to cm


HGCalTBTrack::HGCalTBTrack(void)
{
}

HGCalTBTrack::HGCalTBTrack(const float *v):
	_vertex(MicronToCM * v[1], MicronToCM * v[2], MicronToCM * 0.),
	_momentum(v[3], v[4], 1.)
{
//	assert(size == getSize()); ///\todo to be turned into exception

	chi2	= v[0];
	m_x_err = v[5];
	m_y_err = v[6];

	return;

}

std::vector<float> HGCalTBTrack::getRaw(void)
{
	std::vector<float> v;

	v.push_back(chi2);
	v.push_back(_vertex.X());
	v.push_back(_vertex.Y());
	v.push_back(_momentum.X());
	v.push_back(_momentum.Y());
	v.push_back(m_x_err);
	v.push_back(m_y_err);

	return v;
}


HGCalTBTrack::Point HGCalTBTrack::pointAt(double z) const
{
	return _vertex + _momentum * z;
}
