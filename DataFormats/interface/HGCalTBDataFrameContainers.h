#ifndef HGCalTBDataFrameContainers_h_included
#define HGCalTBDataFrameContainers_h_included 1

#include "DataFormats/Common/interface/DataFrameContainer.h"
#include "HGCal/DataFormats/interface/SKIROC2DataFrame.h"
#include "HGCal/DataFormats/interface/SKIROCParameters.h"

template <class Digi>
class HGCalDataFrameContainer : public edm::DataFrameContainer
{
public:
	typedef edm::DataFrameContainer::size_type size_type;

	HGCalDataFrameContainer(int nsamples_per_digi = SKIROC::MAXSAMPLES, int isubdet = 0) :
		edm::DataFrameContainer(nsamples_per_digi * Digi::WORDS_PER_SAMPLE + Digi::HEADER_WORDS + Digi::FLAG_WORDS, isubdet) { }
	void swap(DataFrameContainer& other)
	{
		this->DataFrameContainer::swap(other);
	}

	//helpful accessors
	using edm::DataFrameContainer::push_back;
	Digi backDataFrame()
	{
		return Digi(this->back());
	}
	int samples() const
	{
		return int((stride() - Digi::HEADER_WORDS - Digi::FLAG_WORDS) / Digi::WORDS_PER_SAMPLE);
	}
	void addDataFrame(DetId detid, const uint16_t* data)
	{
		push_back(detid.rawId(), data);
	}
	void addDataFrame(DetId detid)
	{
		std::vector<uint16_t> dummy_data(stride(), 0);
		push_back(detid.rawId(), &(dummy_data[0]));
	}
	void push_back(const Digi& digi)
	{
		push_back(digi.id(), digi.begin());
	}
};

typedef HGCalDataFrameContainer<SKIROC2DataFrame> SKIROC2DigiCollection; /// \todo TO BE REMOVED


#endif // HGCalTBDataFrameContainers_h_included
