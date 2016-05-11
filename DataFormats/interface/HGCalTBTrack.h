#ifndef HGCALTBTRACK_H
#define HGCALTBTRACK_H

#include <vector>
#include <cstddef>

#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"


/** \class HGCal/DataFormats/interface/HGCalTBTrack.h HGCalTBTrack.h HGCalTBTrack
 * \brief Track from FNAL TB Telescope
 *
 */

class HGCalTBTrack
{
public:
	HGCalTBTrack(void); ///< default constructor
	HGCalTBTrack(const float *raw); ///< constructor from FEDRawData: size is the number of floats in the array
	std::vector<float> getRaw(void); ///< returns a vector of floats to be saved into FEDRawData

	static inline unsigned int getSize(void)
	{
		return 7;
	};
	static inline size_t getSizeof(void)
	{
		return sizeof(float) * getSize() ;
	};

	/// spatial vector
	typedef math::XYZVector Vector;

	/// point in the space
	typedef math::XYZPoint Point;

	inline const Point& vertex() const
	{
		return _vertex;
	};
	inline const Vector& momentum() const
	{
		return _momentum;
	};

	Point pointAt(double z) const;
private:
	Point _vertex; ///< intercept on the telescope plane (z = 0)
	Vector _momentum; ///< m_x(z), m_y(z), 1

	float chi2; ///< chi2 of the track
	float m_x_err; ///< uncertainty on x-direction
	float m_y_err; ///< uncertainty on y-direction
};

#endif

