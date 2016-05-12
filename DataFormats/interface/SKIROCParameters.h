#ifndef SKIROC_PARAMETERS_H
#define SKIROC_PARAMETERS_H

namespace SKIROC
{
static const unsigned char NCHANNELS = 64; ///< number of channels read by one SKIROC
static const unsigned char MAXSAMPLES = 2; ///< one sample for high gain and one for low gain // 15 is the real max of the skiroc
static const unsigned int NLAYERS  = 1; ///< number of layers \todo this should not be hard coded here
}


#endif

