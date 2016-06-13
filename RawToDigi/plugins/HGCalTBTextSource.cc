#include <iostream>
#include "HGCal/RawToDigi/plugins/HGCalTBTextSource.h"
using namespace std;
int runcounter = 0;
bool HGCalTBTextSource::readLines()
{
        int counter = 0;
	m_lines.clear();
	char buffer[1024];
/*
        if(m_file.peek()!='C'){
              throw cms::Exception("MismatchInputSource") << "#" << m_file.peek() << "#";
          }

        buffer[0] = 0;
        m_file.getline(buffer, 1000);
        if( sscanf(buffer, " RUN: %u", &m_run_tmp) != 1) return false;
*/

        while (!feof(m_file) && runcounter < 600) {
//	while (!feof(m_file) && runcounter < 1001) {
		buffer[0] = 0;
		fgets(buffer, 1000, m_file);
//                if (strstr(buffer, "CHIP")) counter++;
//                if (counter == 2) break; // done with this event(2 SKIROCS)
//		if (strstr(buffer, "DONE")) break; // done with this event!
//		if (buffer[0] != '0' && buffer[1] != 'x') continue;
                //if (buffer[0] != ' ') continue;
                
                if (strstr(buffer, "STARTING")) continue;
                if (strstr(buffer, "Board")) continue;
                if(strstr(buffer, "Event")) continue;
                counter++;
		if(counter <= 64) m_lines.push_back(buffer);
                if(counter == 68){
//                   if((runcounter%150) == 0) m_event = 1;
//                   else m_event++;  
                   runcounter++;
                   break; 
                  } 
	}
	return !m_lines.empty();
}

void HGCalTBTextSource::produce(edm::Event & e)
{
	std::auto_ptr<FEDRawDataCollection> bare_product(new  FEDRawDataCollection());


	// here we parse the data
	std::vector<uint16_t> skiwords;
	// make sure there are an even number of 32-bit-words (a round number of 64 bit words...
	if (m_lines.size() % 2) {
		skiwords.push_back(0);
		skiwords.push_back(0);
	}
	for (std::vector<std::string>::const_iterator i = m_lines.begin(); i != m_lines.end(); i++) {
		uint32_t a, b, c, d;
		sscanf(i->c_str(), "%x %x %x %x", &a, &b, &c, &d);
//                cout<<endl<<dec<<a<<dec<<" "<<b<<hex<<" "<<c<<hex<<" "<<d<<endl;
		skiwords.push_back(uint16_t(c >> 16));
                skiwords.push_back(uint16_t(d >> 16));
		skiwords.push_back(uint16_t(c));
                skiwords.push_back(uint16_t(d));
	}

	FEDRawData& fed = bare_product->FEDData(m_sourceId);
	size_t len = sizeof(uint16_t) * skiwords.size();
//        cout<<endl<<" Len= "<<len<<endl;
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
