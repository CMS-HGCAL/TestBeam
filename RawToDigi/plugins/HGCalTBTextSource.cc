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
	m_lines.clear();
	char buffer[1024];
	unsigned int triggerID=0;

	//unsigned int length, bcid;
#ifdef DEBUG
	std::cout << "[DEBUG] Readline" << std::endl;
#endif

	if(m_file.peek()!='C'){
		//cms::LogError("InputSource") << "Input file format does not match: reading character #" << m_file.peek() << "#";
		throw cms::Exception("MismatchInputSource") << "#" << m_file.peek() << "#";
	}

	// read the first line
	buffer[0] = 0;
	m_file.getline(buffer, 1000);
	if( sscanf(buffer, "CHIP %u TRIG: %x TIME: %x RUN: %u EV: %u", &m_sourceId, &triggerID, &m_time, &m_run, &m_event) != 5) return false;
	++m_event;
	//std::cout << triggerID << "\t" << (triggerID & 0xF0000000) << std::endl;

	m_time = m_time >> 2;
#ifdef DEBUGTIME
	std::cout << m_run << "\t" << m_event << "\t" << m_time << "\n";
#endif

	assert( (triggerID & 0xF0000000) == 0x80000000 ); // check if the skiroc is fine 
	
	while ( m_file.peek()!='C' && m_file.good()) {
	  buffer[0] = 0;
	  m_file.getline(buffer, 1000);
//	  assert(buffer[1]=='x');
#ifdef DEBUG
	  std::cout << m_sourceId << "\t" << m_event << "\t" << m_time << "\t" << buffer << "\n"; // buffer has a \n
#endif	  
	  m_lines.push_back(buffer);
	}

	return !m_lines.empty();
}

void HGCalTBTextSource::produce(edm::Event & e)
{
	std::auto_ptr<FEDRawDataCollection> bare_product(new  FEDRawDataCollection());


	// here we parse the data
	std::vector<uint16_t> skiwords;

	for (std::vector<std::string>::const_iterator i = m_lines.begin(); i != m_lines.end(); i++) {
		uint32_t a, b;
		sscanf(i->c_str(), "%x %x", &a, &b);
		if(a==0) continue; // skip the first line because it is the trigger number
		skiwords.push_back(uint16_t(b >> 16));
		skiwords.push_back(uint16_t(b));
	}

	FEDRawData& fed = bare_product->FEDData(m_sourceId);
	size_t len = sizeof(uint16_t) * skiwords.size();
	fed.resize(len);
	memcpy(fed.data(), &(skiwords[0]), len);

	e.put(bare_product);
}



void HGCalTBTextSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setComment("TEST");
	desc.addUntracked<int>("run", 101);
	desc.addUntracked<std::vector<std::string> >("fileNames");
	descriptions.add("source", desc);
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBTextSource);
