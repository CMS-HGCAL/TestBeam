#include <iostream>
#include "HGCal/RawToDigi/plugins/HGCalTBTextSource.h"
using namespace std;

bool HGCalTBTextSource::readLines()
{
        int counter = 0;
	m_lines.clear();
	char buffer[1024];
	while (!feof(m_file)) {
		buffer[0] = 0;
		fgets(buffer, 1000, m_file);
                counter++; 
//                if (strstr(buffer, "CHIP")) counter++;
//                if (counter == 2) break; // done with this event(2 SKIROCS)
//		if (strstr(buffer, "DONE")) break; // done with this event!
//		if (buffer[0] != '0' && buffer[1] != 'x') continue;
                if (buffer[0] != ' ') continue;
                if (strstr(buffer, "  0  0x")) continue;
		m_lines.push_back(buffer);
                if(counter == 132) break; 
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
