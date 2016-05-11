#include "HGCal/CondObjects/interface/HGCalTBNumberingScheme.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

static const int CELLSLOTS_PER_SENSOR = 133 + 2;

HGCalTBNumberingScheme HGCalTBNumberingScheme::the_scheme;

size_t HGCalTBNumberingScheme::rangeFor(uint64_t scheme) const
{
	if (scheme == 0) return 28 * CELLSLOTS_PER_SENSOR;
	else return 0;
}

static const uint8_t numbering_128[] = {
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   0,   1, 255, 255, 255,
	255, 255, 255, 255, 255, 255, 255, 255,   2,   3,   4,   5,   6, 255, 255,
	255, 255, 255, 255, 255, 255,   7,   8,   9,  10,  11,  12,  13,  14, 255,
	255, 255, 255, 255,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
	255, 255, 255,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,
	255, 255, 255,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48, 255,
	255, 255,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60, 255,
	255, 255,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71, 255, 255,
	255,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83, 255, 255,
	255,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94, 255, 255, 255,
	95,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 255, 255, 255,
	107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 255, 255, 255, 255,
	255, 118, 119, 120, 121, 122, 123, 124, 125, 255, 255, 255, 255, 255, 255,
	255, 255, 126, 127, 128, 129, 130, 255, 255, 255, 255, 255, 255, 255, 255,
	255, 255, 255, 131, 132, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255
};

/// \todo need documentation!!!
size_t HGCalTBNumberingScheme::denseIndexFor(uint32_t rawDetId, uint64_t scheme) const
{
	HGCalTBDetId id(rawDetId);
	if (scheme != 0) return HGCalCondObjectNumberingScheme::INVALID;

	size_t idx = (id.layer() - 1) * CELLSLOTS_PER_SENSOR;
	int linear = (id.iu() + 7) * 15 + id.iv() + 7;
	if (id.cellType() == HGCalTBDetId::kCellTypeCalibInner) {
		if (id.iv() > 0) idx += 133;
		else idx += 134;
	} else idx += numbering_128[linear];

	return idx;
}
