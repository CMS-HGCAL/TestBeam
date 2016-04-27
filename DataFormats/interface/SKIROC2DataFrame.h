#ifndef SKIROC2DATAFRAME_H_INCLUDED
#define  SKIROC2DATAFRAME_H_INCLUDED 1

#include <stdint.h>
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "DataFormats/Common/interface/DataFrame.h"

class SKIROC2DataFrame
{
public:

	static const int WORDS_PER_SAMPLE = 2;
	static const int HEADER_WORDS = 1;
	static const int FLAG_WORDS = 1;

	SKIROC2DataFrame() { }
	SKIROC2DataFrame(edm::DataFrame const & df) : m_data(df) { }

	class Sample
	{
	public:
		Sample(const edm::DataFrame& frame, edm::DataFrame::size_type i) : frame_(frame), i_(i) { }
		static const int MASK_ADCTDC = 0x0FFF;
		static const int MASK_HIT = 0x1000;
		static const int MASK_GAIN = 0x2000;

		enum {
			ADCLOW_SHIFT = 0,
			ADCHIGH_SHIFT,
			TDC_SHIFT
		};

		int adcLow() const // adc values in 12 bits, so no reason to use unsigned types
		{
			return frame_[i_ + ADCLOW_SHIFT] & MASK_ADCTDC;
		}
		int adcHigh() const
		{
			return frame_[i_ + ADCHIGH_SHIFT] & MASK_ADCTDC;
		}

		int tdc() const
		{
			return frame_[i_ + TDC_SHIFT] & MASK_ADCTDC;
		}
	private:
		const edm::DataFrame& frame_;
		edm::DataFrame::size_type i_;
	};

	void copyContent(const SKIROC2DataFrame& src);

	/// Get the detector id
	HGCalTBDetId detid() const
	{
		return HGCalTBDetId(m_data.id());
	}

	edm::DataFrame::id_type id() const
	{
		return m_data.id();
	}
	/// more accessors
	edm::DataFrame::size_type size() const
	{
		return m_data.size();
	}
	/// iterators
	edm::DataFrame::iterator begin()
	{
		return m_data.begin();
	}
	edm::DataFrame::iterator end()
	{
		return m_data.end();
	}
	edm::DataFrame::const_iterator begin() const
	{
		return m_data.begin();
	}
	edm::DataFrame::const_iterator end() const
	{
		return m_data.end();
	}

	/// total number of samples in the digi
	int samples() const
	{
		return (size() - HEADER_WORDS - FLAG_WORDS) / WORDS_PER_SAMPLE;
	}
	/// get the sample
	inline Sample operator[](edm::DataFrame::size_type i) const
	{
		return Sample(m_data, i * WORDS_PER_SAMPLE + HEADER_WORDS);
	}
	/// set the sample contents
	void setSample(edm::DataFrame::size_type isample, int adcLow, int adcHigh, int tdc);
	/// get the flag word
	uint16_t flags() const
	{
		return m_data[size() - 1];
	}
	/// set the flag word
	void setFlags(uint16_t v);

private:
	edm::DataFrame m_data;
};

std::ostream& operator<<(std::ostream&, const SKIROC2DataFrame&);

#endif // SKIROC2DATAFRAME_H_INCLUDED
