#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <iostream>
#include <stdio.h>
#include <fstream>

/**
   \page TXTFORMAT_PAGE Data txt input format:
   \tableofcontents
   Return to the main page: \ref index
   \section TXTFORMAT Original data txt format:
    - Spill header containing spill time, run number and event number in the spill
	\verbatim
	STARTING SPILL READ AT TIME (1us): 0x55DC7602 RUN: 802 EVENT: 150
	\endverbatim
	- Board header with the FMC-IO identification number
	\verbatim
	Board header: on FMC-IO 11, trig_count in mem= 150, sk_status = 1
	\endverbatim
	- Event header with event number, timestamp of the event and the global trigger time
	\verbatim
	Event header for event 0 with (200ns) timestamp 0xB803BA11  global tts(us)  0x00000000 and CKOV= 0
	\endverbatim
	- 68 lines with: eventNumber, channelID, 32bit word formed by two 16bit words for low gain values of two channels, 32bit word formed by two 16bit words for high gain vaues of two channels
	the 68 lines have then the information from two skirocs
	\verbatim
	0 0 0x11B411B8   0x11941188
	\endverbatim

	The original txt format is not suitable for being processed by CMSSW, events are not ordered. For each board all the events are dumped.

	\section NEWTXTFORMAT_ Rearranged data txt format:
	In order to be able to process the data with CMSSW, we need to rearrange the txt in order to have for each event all the boards.

	This is done with the script \verbatim ./scripts/rearrangeTxtFile.sh <originalFile.txt> <output_directory> \endverbatim

	This step is done centrally and files are available in ...

	The txt file is processed by HGCalTBTextSource

	\section EXAMPLE Usage example in your python cfg:
	\code
	process.source = cms.Source("HGCalTBTextSource",
	  fileNames=cms.untracked.vstring("file:myfile1.txt"),
	)
	\endcode

	Return to the main page: \ref index

	\example test_cfg.py

*/

/**
 * \class HGCalTBTextSource HGCal/RawToDigi/plugins/HGCalTBTextSource.h
 *
 * \brief Convert data from txt file to FEDRawData
 *
 * \details
 * For info about the txt input format see: \ref TXTFORMAT_PAGE
 * \author Shervin Nourbakhsh (UMN)
 * \author Rajdeep Mohan Chatterjee (UMN)
 * \author Jeremy Mans (UMN)
 *
 * \todo fedID is taken from cfg, it should be taken from the txt file
 *
 * DATA in the FEDRAWDATA in the following order:
 * BOARD, SKIROC, CHANNEL
 * so [Board0 [Skiroc0 [Channel0][Channel1][Channel2]...][Skiroc1 [Channel0][Channel1][Channel2]...]]
 */
class HGCalTBTextSource : public edm::ProducerSourceFromFiles
{

public:
	explicit HGCalTBTextSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
		m_file(0)
	{
		m_sourceId = pset.getUntrackedParameter<int>("fed", 1000); /// \todo read from file
		produces<FEDRawDataCollection>();
	}

	virtual ~HGCalTBTextSource()
	{
		if (m_file != 0) fclose(m_file);
		m_file = 0;
	}

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);

	/** DATA in the FEDRAWDATA in the following order:
	 * BOARD, SKIROC, CHANNEL
	 * so [Board0 [Skiroc0 [Channel0][Channel1][Channel2]...][Skiroc1 [Channel0][Channel1][Channel2]...]]
	 */
	virtual void produce(edm::Event & e); ///test
	bool readHeader(void);
	bool readLines(void);

	std::array< std::array< std::vector < unsigned int> ,  MAXSKIROCS_PER_BOARD >, MAXLAYERS > m_lines;
	FILE* m_file;
	unsigned int m_time;
	unsigned int m_event, m_run, m_spill, max_boards;
	int m_sourceId;
};
