#include <iostream>
#include "HGCal/RawToDigi/plugins/HGCalTBTextSource.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"
using namespace std;
unsigned int runcounter = 0;
int runflag = 0;
int DangerFlag = 0;// Catch Danger and push back -10000
int trigcountperspillflag = 0;
unsigned int Events_Per_Spill = 0;
unsigned int spillcounter = 0;
unsigned int tmp_event = 0;
//unsigned int Number_Of_Events_Per_Spill = 150;
//unsigned int Number_Of_Spills = 8;
int Number_Of_SKIROC_Data_Words = 64;
int Number_Of_SKIROC_Words = 68;

string buffer1 = "0 0 0xFFFFFFFF   0xFFFFFFFF";
bool HGCalTBTextSource::readLines()
{
	m_lines.clear();
	char buffer[1024];
	int counter = 0;
	if(runcounter == 0 && runflag == 0) {
		buffer[0] = 0;
		fgets(buffer, 1000, m_file);
		if(sscanf(buffer, "STARTING SPILL READ AT TIME (1us): %x RUN: %u", &m_time_tmp, &m_run_tmp) != 2 ) return false;
		m_run = m_run_tmp;
		runflag = 1;
	}
	if(runcounter == EVENTSPERSPILL) {
		runcounter = 0;
		trigcountperspillflag++;
	}

	if(trigcountperspillflag == MAXLAYERS) {
		spillcounter++;
		trigcountperspillflag = 0;
	}

	/*
	         if(runcounter == 0 && trigcountperspillflag == 0){
	               buffer[0] = 0;
	               fgets(buffer, 1000, m_file);
	               if(sscanf(buffer,"Board header: on FMC-IO 2, trig_count in mem= %u, sk_status = 1",&Events_Per_Spill) !=1 ) return false;
	               trigcountperspillflag = 1;
	               if(Events_Per_Spill < 140) return false;
	             }
	*/

//        if( sscanf(buffer, " RUN: %u", &m_run) != 1) return false;
//        sscanf(buffer, " RUN:%u", &m_run);
//        cout<<endl<<m_run<<endl;
	while (!feof(m_file) && spillcounter < m_nSpills) {
//	while (!feof(m_file) && runcounter < 1001) {
		buffer[0] = 0;
		fgets(buffer, 1000, m_file);
		if (strstr(buffer, "STARTING")) continue;
		if (strstr(buffer, "Board")) continue;
//                if(strstr(buffer, "Event")) continue;


		/*
		                if(strstr(buffer,"CKOV= 1")){
		                    for(int iii=1; iii<=64;iii++){
		 			m_lines.push_back(buffer1);
		                        if(iii==64) counter=67;
		                       }
		                   }

		                if(strstr(buffer,"CKOV= 1")) cout<<endl<<"CKOV= 1"<<endl;
		                else cout<<endl<<"CKOV= 0"<<endl;
		*/

		if(sscanf(buffer, "Event header for event %x with (200ns) timestamp %x", &tmp_event, &m_time) == 2) {
			DangerFlag = 0;
			continue;
		}

		if(sscanf(buffer, "DANGER: stale event warning! Event header for event %x with (200ns) timestamp %x", &tmp_event, &m_time) == 2){
			DangerFlag = 1;
			continue;
		}	

//                if(strstr(buffer, "Event")) continue;
		counter++;
		if((DangerFlag == 0) && runcounter < EVENTSPERSPILL && counter <= Number_Of_SKIROC_Data_Words) m_lines.push_back(buffer);
                if((DangerFlag == 1) && runcounter < EVENTSPERSPILL && counter <= Number_Of_SKIROC_Data_Words) m_lines.push_back(buffer1);
		if(runcounter < EVENTSPERSPILL && counter == Number_Of_SKIROC_Words) {
			runcounter++;
//                   cout<<endl<<" Reached Here 1"<<endl;
			break;
		}

	}

	return !m_lines.empty();
}

void HGCalTBTextSource::produce(edm::Event & e)
{
	std::auto_ptr<FEDRawDataCollection> bare_product(new  FEDRawDataCollection());

//        cout<<endl<<"Size= "<<m_lines.size()<<endl;
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
//                cout<<endl<<" Here "<<dec<<a<<dec<<" "<<b<<hex<<" "<<c<<hex<<" "<<d<<endl;
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
	desc.addUntracked<unsigned int>("nSpills", 6);
	descriptions.add("source", desc);
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBTextSource);
