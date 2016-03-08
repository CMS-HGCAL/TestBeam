#include "HGCal/CondObjects/interface/HGCalCondObjectContainer.h"

size_t HGCalCondObjectContainerBase::indexOf(DetId id) const
{
	if (!p_numbering) return HGCalCondObjectNumberingScheme::INVALID;
	else return p_numbering->denseIndexFor(id.rawId(), m_numbering_code);
}
