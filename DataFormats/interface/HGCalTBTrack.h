#ifndef HGCALTBTRACK_H
#define HGCALTBTRACK_H

#include <vector>
#include <cstddef>

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

	static inline unsigned int getSize(void) { return 7;};
	static inline size_t getSizeof(void) { return sizeof(float) * getSize() ;};
private:
	float chi2; ///< chi2 of the track
	float x0;   ///< intercept: x on the telescope plane
	float y0;   ///< intercept: y on the telescope plane
	float m_x;  ///< x direction from telescope hits
	float m_y;  ///< y direction from telescope hits
	float m_x_err; ///< uncertainty on x-direction
	float m_y_err; ///< uncertainty on y-direction

};

#endif

