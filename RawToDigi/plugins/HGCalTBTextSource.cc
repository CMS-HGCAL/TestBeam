#include <iostream>
#include <sstream>
#include <string>
#include "stdlib.h"
#include "HGCal/RawToDigi/plugins/HGCalTBTextSource.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#define TXTLINES 69

//#define DEBUG
/**
RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=09	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=12	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=14	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=17	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=18	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=19	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=22	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=23
T 0x0072D911 0xBEC20B15	T 0x0072D93F 0xBED19F54	T 0x0072D96E 0xBED19F56	T 0x0072D99D 0x3A11F91F	T 0x0072D9CC 0x40DACEB1	T 0x0072D9FA 0xC5EAE85C	T 0x00000000 0xBED19F52	T 0x00000000 0x57B86176
0x118F11BD 0x1081108D	0x11A511B7 0x118A108C	0x119311F2 0x119F1195	0x109E109E 0x1083108D	0x118D11F2 0x11FD1198	0x11881195 0x118D1196	0x11911190 0x1193118A	0x1088109F 0x118C118A
*/

bool HGCalTBTextSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType& evType)
{

	if (!(fileNames().size())) return false; // need a file...
	if (m_file == 0) {
		std::string name = fileNames()[0].c_str(); /// \todo FIX in order to take several files
		if (name.find("file:") == 0) name = name.substr(5);
		m_file = fopen(name.c_str(), "r");
		if (m_file == 0) {
			throw cms::Exception("FileNotFound") << "Unable to open file " << name;
		}
		m_event = 0;
	}
	if (feof(m_file)) return false;
	// read the header
	while(readHeader() == false) {
		if (feof(m_file)) return false;
		// have to skip the data lines and return false
		readLines();
	}

	// now try to read from the file
	if (!readLines()) return false;

	id = edm::EventID(m_run, m_spill, m_event);

	time = (edm::TimeValue_t) m_time;
	return true;
}

// sets m_run, m_spill, m_even, m_time
///\todo read the number of boards
bool HGCalTBTextSource::readHeader()
{
	char buffer[1024]; // typically 605 chars
	buffer[0] = 0;     //init buffer
	fgets(buffer, 1000, m_file); // read the header line
	std::string b = buffer;
#ifdef DEBUG
	std::cout << "------------------------------\n";
	std::cout << b << std::endl;
#endif
	if(b.find("DANGER=true") != std::string::npos) return false;
//RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=09	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=12	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=14	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=17	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=18	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=19	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=22	RUN=000880	SPILL=01	EVENT=000000	GLOBALTIME=0x01B6EAF2	BOARD=23

	if(sscanf(buffer, "RUN=%u\tSPILL=%u\tEVENT=%u\tGLOBALTIME=%x", &m_run, &m_spill, &m_event, &m_time) != 4) {
		return false; // what to do?
	}

	// here count the number of boards for this event

	return true;
}


// each time it returns true, then the event is complete and can be saved in EDM format
bool HGCalTBTextSource::readLines()
{
	max_boards = 0;
	for( auto& board : m_lines) {
		for( auto& skiroc : board){
			skiroc.clear();
		}
	}

	char buff[1024], buff_SK[1024];
	buff[0] = 0;
	unsigned int data_sk;

	// loop over all the lines
	for(unsigned int i = 0; i < TXTLINES && !feof(m_file); ++i) {
		buff[0] = 0;
		fgets(buff, 1000, m_file);
		if((i < 1) || (i > 64)) continue;// Only these are data words, the rest dont interest us for now
		// loop over one line of the text file
		std::string b = buff;
		std::istringstream buffer(b);
		unsigned int board_counter = 0;

		while( buffer.peek() != '\n' && buffer.good()) {
			//continue for all the boards
			for(size_t iSkiroc = 0; iSkiroc < MAXSKIROCS_PER_BOARD; ++iSkiroc) {
				// read the data of the one skiroc of one board
				buffer >> buff_SK;
				data_sk = strtoul(buff_SK, NULL, 0);
				m_lines[board_counter][iSkiroc].push_back(data_sk);
			}
			++board_counter;
			if(board_counter > max_boards) max_boards = board_counter;
		}

	}
	return true;
}

/** DATA in the FEDRAWDATA in the following order:
 * BOARD, SKIROC, CHANNEL
 * so [Board0 [Skiroc0 [Channel0][Channel1][Channel2]...][Skiroc1 [Channel0][Channel1][Channel2]...]]
 */
void HGCalTBTextSource::produce(edm::Event & e)
{
	std::auto_ptr<FEDRawDataCollection> bare_product(new  FEDRawDataCollection());
	// here we parse the data
	std::vector<uint16_t> skiwords;
	// make sure there are an even number of 32-bit-words (a round number of 64 bit words...

	for (size_t i_board = 0 ; i_board < max_boards; ++i_board) {
		for(size_t i_skiroc = 0; i_skiroc < MAXSKIROCS_PER_BOARD; ++i_skiroc) {
			auto& board = m_lines[i_board][i_skiroc];
			for (auto& skiword : board) {
				skiwords.push_back(skiword >> 16);
				skiwords.push_back(skiword);
			}
		}
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
