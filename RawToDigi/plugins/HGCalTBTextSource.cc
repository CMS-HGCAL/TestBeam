#include "HGCal/RawToDigi/plugins/HGCalTBTextSource.h"
//#define DEBUG
//#define DEBUGTIME

#ifdef DEBUG
#include <iostream>
#endif

#ifdef DEBUGTIME
#include <iostream>
#endif

void HGCalTBTextSource::openFile(edm::FromFiles& files, std::ifstream& file)
{

	if(file.is_open() && !file.eof()) return; // continue to use the same file if the end is not reached
	if(file.is_open()){ //the end of file is reached
		file.close();
		files.incrementFileIndex();
	}

	std::string name = files.fileNames()[files.fileIndex()];
	
	if (name.find("file:") == 0) name = name.substr(5);
#ifdef DEBUG
	std::cout << "[DEBUG] Opening file: " << name << std::endl;
#endif
	
	file.open(name);
	if (file.fail()) {
		throw cms::Exception("FileNotFound") << "Unable to open file " << name;
	}
}


bool HGCalTBTextSource::readLines()
{
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
			_triggerID = triggerID_tmp;

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

bool HGCalTBTextSource::readTelescopeLines()
{
	if(!_telescope_words.empty()) return true; // data have not been used, don't read more
	if(_telescopeFile.peek() && !_telescopeFile.good()) return false;

	char buffer[1024];
	buffer[0] = 0;
	_telescopeFile.getline(buffer, 1000);

	parseAddTelescopeWords(_telescope_words, buffer);

	return !_telescope_words.empty();
}


bool HGCalTBTextSource::setRunAndEventInfo(edm::EventID& id, edm::TimeValue_t& time, edm::EventAuxiliary::ExperimentType&)
{
	
	openFile(*this, _hgcalFile); // 
	openFile(_telescopeFiles, _telescopeFile);

	if(!readLines()) return false; // readlines is here because the run and event info are contained in the file 
	if(!readTelescopeLines()) return false;
	
	id = edm::EventID(m_run, 1, m_event);
	
	// time is a hack
	time = (edm::TimeValue_t) m_time;
#ifdef DEBUGTIME
	std::cout << m_time << "\t" << time << std::endl;
#endif
	return true;
}

void HGCalTBTextSource::produce(edm::Event & e)
{
	std::auto_ptr<FEDRawDataCollection> bare_product(new  FEDRawDataCollection());

	FEDRawData& fed = bare_product->FEDData(m_sourceId);
	size_t len = sizeof(uint16_t) * m_skiwords.size();
	fed.resize(len);
	memcpy(fed.data(), &(m_skiwords[0]), len);

	fed = bare_product->FEDData(_TELESCOPE_FED_ID_);
	if(_triggerID==t_triggerID){ // empty FED if no data are available for the triggerID
		len = sizeof(float) * telescope_words.size();
		fed.resize(len);
		memcpy(fed.data(), &(telescope_words[0]), len);
	}

	// words are reset, only if the vectors are empty new lines are going to be read
	m_skiwords.clear();
	_telescope_words.clear();

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

void HGCalTBTextSource::parseAddTelescopeWords(std::vector<float>& telescope_words, std::string i){
//# UTC_Time_Stamp,Trigger Number,Number_Of_Tracks,Chi2/NDF,X_Intercept,Y_Intercept,X_Slope,Y_Slope,X_Slope_Error,Y_Slope_Error
//1459822031,7,1,1.0206740673,19609.2010917080,13026.1630167309,-0.0001089257,-0.0002085995,0.0000209990,0.0000139143
	unsigned int time, nTracks;
	sscanf(i.c_str(), "%u,%u,%u,%f,%f,%f,%f,%f,%f", &time, &t_triggerID, &nTracks, &chi2, &x0, &y0, &m_x, &m_y, &m_x_err, &m_y_err);
	telescope_words.push_back(chi2);
   	telescope_words.push_back(x0);
	telescope_words.push_back(y0);
	telescope_words.push_back(m_x);
	telescope_words.push_back(m_y);
	telescope_words.push_back(m_x_err);
	telescope_words.push_back(m_y_err);

}

void HGCalTBTextSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setComment("TEST");
	desc.addUntracked<std::vector<std::string> >("fileNames");
	desc.addUntracked<std::vector<std::string> >("telescopeFiles");
	descriptions.add("source", desc);
}


#include "FWCore/Framework/interface/InputSourceMacros.h"

DEFINE_FWK_INPUT_SOURCE(HGCalTBTextSource);
