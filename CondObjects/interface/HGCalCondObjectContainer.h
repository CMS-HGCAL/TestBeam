#ifndef HGCAL_CONDOBJECTS_HGCALCONDOBJECTCONTAINER_H
#define HGCAL_CONDOBJECTS_HGCALCONDOBJECTCONTAINER_H 1

#include <vector>
#include "DataFormats/DetId/interface/DetId.h"


class HGCalCondObjectNumberingScheme
{
public:
	static const size_t INVALID = (size_t) - 1;
	virtual size_t rangeFor(uint64_t scheme) const = 0;
	virtual size_t denseIndexFor(uint32_t rawDetId, uint64_t scheme) const = 0;
	virtual size_t denseIndexFor(DetId rawDetId, uint64_t scheme) const
	{
		return denseIndexFor(rawDetId.rawId(), scheme);
	}
};


/** \class HGCalCondObjectContainerBase
  *
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
class HGCalCondObjectContainerBase
{
public:
	HGCalCondObjectContainerBase(const HGCalCondObjectNumberingScheme* scheme, uint64_t ischeme) : m_numbering_code(ischeme), p_numbering(scheme) { }
	HGCalCondObjectContainerBase() : m_numbering_code(0), p_numbering(0) { }

	void setNumberingScheme(const HGCalCondObjectNumberingScheme* scheme)
	{
		p_numbering = scheme;
	}
	void setNumberingScheme(const HGCalCondObjectNumberingScheme* scheme, uint64_t code)
	{
		p_numbering = scheme;
		m_numbering_code = code;
	}
	const HGCalCondObjectNumberingScheme* getNumberingScheme() const
	{
		return p_numbering;
	}
	uint64_t schemeCode() const
	{
		return m_numbering_code;
	}

	bool exists(DetId id) const
	{
		return indexOf(id) != HGCalCondObjectNumberingScheme::INVALID;
	}
protected:
	size_t indexOf(DetId id) const;
private:
	uint64_t m_numbering_code;
	const HGCalCondObjectNumberingScheme* p_numbering; // volatile
};

/** \class HGCalCondObjectContainer
  *
  * $Date: $
  * $Revision: $
  * \author J. Mans - Minnesota
  */
template <class Payload>
class HGCalCondObjectContainer : public HGCalCondObjectContainerBase
{
public:
	struct Item {
		DetId id;
		Payload value;
	};

	HGCalCondObjectContainer(const HGCalCondObjectNumberingScheme* scheme, uint64_t ischeme) : HGCalCondObjectContainerBase(scheme, ischeme), m_payload(scheme->rangeFor(ischeme))
	{
	}
	HGCalCondObjectContainer() : HGCalCondObjectContainerBase() { }

	const Item* get(DetId id) const
	{
		size_t i = indexOf(id);
		return (i == HGCalCondObjectNumberingScheme::INVALID) ? (0) : (&(m_payload[i]));
	}
	size_t size() const
	{
		return m_payload.size();
	}
	const Item& get(size_t i) const
	{
		return m_payload[i];
	}
	void set(DetId id, const Payload& value)
	{
		size_t i = indexOf(id);
		if (i != HGCalCondObjectNumberingScheme::INVALID) {
			m_payload[i].id = id;
			m_payload[i].value = value;
		}
	}
private:
	std::vector<Item> m_payload;
};

#endif
