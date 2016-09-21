#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Sources/interface/ProducerSourceFromFiles.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <iostream>
#include <stdio.h>
#include <fstream>

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
		m_run(pset.getUntrackedParameter<int>("run", 101)), /// \todo check and read from file?
		m_nSpills(pset.getUntrackedParameter<unsigned int>("nSpills", 15))
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
	bool setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&);
	virtual void produce(edm::Event & e);
	bool readLines(void);

	std::vector<std::string> m_lines;
	FILE* m_file;
	unsigned int m_time, m_time_tmp, m_run_tmp;
	int m_event, m_run;
	unsigned int m_nSpills;
	int m_sourceId;
};
