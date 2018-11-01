#ifndef HGCAL_TB_GEOMETRY_PARAMETERS_H
#define HGCAL_TB_GEOMETRY_PARAMETERS_H

//old stuff to be removed
#define MAXSKIROCS 40
#define MAXLAYERS 6

namespace HGCAL_TB_GEOMETRY
{
  static const uint16_t NUMBER_OF_LAYERS = 6;
  static const uint16_t NUMBER_OF_HEXABOARD = 94;
  static const uint16_t N_SKIROC_PER_HEXA = 4;
  static const uint16_t N_CHANNELS_PER_SKIROC = 64;

  static const int MAXVERTICES = 6;
  static const double DELTA = 0.00001;//Add/subtract delta = 0.00001 to x,y of a cell centre so the TH2Poly::Fill doesnt have a problem at the edges where the centre of a half-hex cell passes through the sennsor boundary line.
}
#endif

