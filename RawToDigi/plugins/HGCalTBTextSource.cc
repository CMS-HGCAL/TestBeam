#include "HGCal/RawToDigi/plugins/HGCalTBTextSource.h"
//#define DEBUG
#define DEBUGTIME

#ifdef DEBUG
#include <iostream>
#endif

#ifdef DEBUGTIME
#include <iostream>
#endif


bool HGCalTBTextSource::readLines()
{
	//m_lines.clear();
	m_skiwords.clear();
	char buffer[1024];
	unsigned int triggerID_tmp=0;

	//unsigned int length, bcid;
#ifdef DEBUG
	std::cout << "[DEBUG] Readline" << std::endl;
#endif

	if(m_file.peek()!='C'){
		//cms::LogError("InputSource") << "Input file format does not match: reading character #" << m_file.peek() << "#";
		throw cms::Exception("MismatchInputSource") << "#" << m_file.peek() << "#";
	}

	for(unsigned int iSkiroc =0 ; iSkiroc < MAXSKIROCS; ++iSkiroc){

		// read the first line
		buffer[0] = 0;
		m_file.getline(buffer, 1000);
		if( sscanf(buffer, "CHIP %u TRIG: %x TIME: %x RUN: %u EV: %u", &m_sourceId_tmp, &triggerID_tmp, &m_time_tmp, &m_run_tmp, &m_event_tmp) != 5) return false;
		if(iSkiroc==0){
			m_event = m_event_tmp;
			m_run = m_run_tmp;
			m_time = m_time_tmp;
			m_sourceId = m_sourceId_tmp;

			++m_event;
			m_time = m_time >> 8;
		}

#ifdef DEBUGTIME
		std::cout << m_run << "\t" << m_event << "\t" << m_time << "\n";
#endif
		
		assert( (triggerID_tmp & 0xF0000000) == 0x80000000 ); // check if the skiroc is fine  ///\todo transform to exception
		
		while ( m_file.peek()!='C' && m_file.good()) {
			buffer[0] = 0;
			m_file.getline(buffer, 1000);
#ifdef DEBUG
			std::cout << m_sourceId << "\t" << m_event << "\t" << m_time << "\t" << buffer << "\n"; // buffer has a \n
#endif	  
			parseAddSkiword(m_skiwords, buffer);
		}

	}
	return !m_skiwords.empty();
}

void HGCalTBTextSource::produce(edm::Event & e)
{
	std::auto_ptr<FEDRawDataCollection> bare_product(new  FEDRawDataCollection());

	FEDRawData& fed = bare_product->FEDData(m_sourceId);
	size_t len = sizeof(uint16_t) * m_skiwords.size();
	fed.resize(len);
	memcpy(fed.data(), &(m_skiwords[0]), len);

	e.put(bare_product);
}


void HGCalTBTextSource::parseAddSkiword(std::vector<uint16_t>& skiwords, std::string i){
	uint32_t a, b;
	sscanf(i.c_str(), "%x %x", &a, &b);
	if(a==0) return;
		
// two words of 16bits for ADC high gain and low gain
	skiwords.push_back(uint16_t(b >> 16));
	skiwords.push_back(uint16_t(b));
}



void HGCalTBTextSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setComment("TEST");
	desc.addUntracked<std::vector<std::string> >("fileNames");
	descriptions.add("source", desc);
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBTextSource);
