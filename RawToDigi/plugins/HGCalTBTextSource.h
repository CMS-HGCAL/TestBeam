#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <stdio.h>

/**
 * \class HGCalTBTextSource HGCal/RawToDigi/plugins/HGCalTBTextSource.h
 *
 * \brief convert data from txt file to FEDRawData
 *
 * \todo replace c-like scanf with c++ versions
 * \todo change run and fed IDs (now are hardcoded)
 */
class HGCalTBTextSource : public edm::ProducerSourceFromFiles
{

public:
	explicit HGCalTBTextSource(const edm::ParameterSet & pset, edm::InputSourceDescription const& desc) :  edm::ProducerSourceFromFiles(pset, desc, true),
		m_file(0),
		m_run(pset.getUntrackedParameter<int>("run", 101)) /// \todo check and read from file?
	{

		m_sourceId = pset.getUntrackedParameter<int>("fed", 1000); /// \todo check and read from file?
		produces<FEDRawDataCollection>();
	}

	virtual ~HGCalTBTextSource()
	{
		if (m_file != 0) fclose(m_file);
		m_file = 0;
	}

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&)
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
		if (!readLines()) return false;

		m_event++; /// \todo get eventID from file?
		id = edm::EventID(m_run, 1, m_event);
		// time is a hack
		edm::TimeValue_t present_time = presentTime(); ///\todo take time from file? how to define the time?
		unsigned long time_between_events = timeBetweenEvents();

		time = present_time + time_between_events;
		return true;
	}
	virtual void produce(edm::Event & e);
	bool readLines();

	std::vector<std::string> m_lines;
	FILE* m_file;
	int m_event, m_run;
	int m_sourceId;
};
