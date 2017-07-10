// -*- C++ -*-
//
// Package:    UserCode/Bad_Spill_Filter
// Class:      Bad_Spill_Filter
// 
/**\class Bad_Spill_Filter Bad_Spill_Filter.cc UserCode/Bad_Spill_Filter/plugins/Bad_Spill_Filter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Fri, 18 Nov 2016 12:57:18 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string>
//
// class declaration
//

using namespace std;

class Bad_Spill_Filter : public edm::stream::EDFilter<> {
   public:
      explicit Bad_Spill_Filter(const edm::ParameterSet&);
      ~Bad_Spill_Filter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      int layers_config_;
      std::string name_CFG_;
      char buffer[1024];
      int m_run, tmp_run, m_spill, tmp_spill;
      int spill_counter = 0;
      int Number_Of_Bad_Runs = 0;	
      const static int Num_BAD_RUNS_CFG = 400;//Choose a number larger than that in CFG1 or 2, the run numbers are read from the bad run/spill file for each configuration and saved in Number_Of_Bad_Runs in the constructor
      int BAD_Runs_CFG[Num_BAD_RUNS_CFG] = {0};

//      int Bad_Runs_CFG1[115] = {1051, 1054, 1055, 1056, 1058, 1061, 1062, 1063, 1064, 1065, 1067, 1068, 1071, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 857, 858, 859, 860, 861, 1026, 1027, 1029, 1031, 1032, 1033, 1034, 1035, 1036, 1038, 1040, 1041, 1043, 1044, 1045, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 937, 938, 939, 940, 941, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 956, 957, 960, 961, 1087, 1088, 1091, 1093, 1094, 1096, 1097, 1098, 1099, 1103, 1104, 1105, 1106, 1107, 865, 866, 867, 868, 869, 872, 873, 875, 876, 877, 878, 880, 881, 882, 883, 888, 889};		
//	int Bad_Runs_CFG2[98] = {1202, 1203, 1204, 1206, 1208, 1214, 1215, 1216, 1217, 1219, 1220, 1221, 1222, 1223, 1179, 1180, 1181, 1185, 1186, 1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 1197, 1199, 1122, 1123, 1124, 1125, 1126, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1137, 1138, 1141, 1144, 1145, 1146, 1150, 1248, 1249, 1250, 1251, 1253, 1254, 1257, 1260, 1261, 1263, 1264, 1267, 1291, 1292, 1294, 1295, 1296, 1298, 1299, 1301, 1302, 1303, 1305, 1307, 1308, 1309, 1227, 1228, 1233, 1235, 1240, 1242, 1243, 1244, 1246, 1247, 1154, 1156, 1157, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1172, 1173};

      std::array< std::vector <int> , Num_BAD_RUNS_CFG> Bad_Run_Spill_Array;

	  FILE* m_file;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Bad_Spill_Filter::Bad_Spill_Filter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
	bool Not_Filled_Once = true;
	int run_counter = 0;
	layers_config_ = iConfig.getParameter<int>("layers_config");

	if(layers_config_ == 1) name_CFG_ = iConfig.getParameter<std::string>("nameCFG1");
	else if(layers_config_ == 2) name_CFG_ = iConfig.getParameter<std::string>("nameCFG2");


	m_file = fopen(name_CFG_.c_str(),"r");
	if (m_file == 0) {
	  if(layers_config_ == 1) cout<<endl << "Unable to open file for cfg1 " << name_CFG_;
	  if(layers_config_ == 2) cout<<endl << "Unable to open file for cfg2" << name_CFG_;
	  exit(0);
	}

	while(!feof(m_file)){

	        fgets(buffer, 1000, m_file);
       		if(sscanf(buffer, "%u  %u", &m_run, &m_spill) != 2) {
	                continue;
	        }

	        if(Not_Filled_Once){
			tmp_run = m_run;
			BAD_Runs_CFG[run_counter] = m_run;
			Bad_Run_Spill_Array[run_counter].push_back(m_spill);
			Not_Filled_Once = false;	
		}

	        if(!Not_Filled_Once){
			if(m_run != tmp_run){
				run_counter++;
				BAD_Runs_CFG[run_counter] = m_run;
				Bad_Run_Spill_Array[run_counter].push_back(m_spill);
		                tmp_run = m_run;
			}
			else{
	                        Bad_Run_Spill_Array[run_counter].push_back(m_spill);
				}
	        }

	        
	}

	Number_Of_Bad_Runs = run_counter;
	cout<<endl<<" Total number of Bad Runs = "<<Number_Of_Bad_Runs<<endl;

}


Bad_Spill_Filter::~Bad_Spill_Filter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
Bad_Spill_Filter::filter(edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;
		
        int runId = event.id().run();
        int spillId = event.luminosityBlock();
	int BadRunFlag = 0;
	int BadSpillFlag = 0;
	int Bad_Run_Location = 0;

	if((layers_config_ != 1) && (layers_config_ != 2) ) return true;

	for(int iii = 0; iii <= Number_Of_Bad_Runs; iii++){
		if(runId == BAD_Runs_CFG[iii]){
			BadRunFlag = 1;
			Bad_Run_Location = iii;
		}
	}

        if(BadRunFlag == 1){
		auto Bad_Run_Spill_List = Bad_Run_Spill_Array[Bad_Run_Location];
		for(auto Bad_Spill_Iterator : Bad_Run_Spill_List){
			if(spillId == Bad_Spill_Iterator){
				if(tmp_spill != spillId) spill_counter = 0;
				if(spill_counter == 0){
					cout<<endl<<" Run, Spill Filtered out : "<<runId<<" , "<<spillId<<endl;
					tmp_spill = spillId;
				}
				spill_counter++;
				BadSpillFlag = 1;
			}
		}
	}

        if(BadSpillFlag == 1) return false;
        return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
Bad_Spill_Filter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
Bad_Spill_Filter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
Bad_Spill_Filter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
Bad_Spill_Filter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Bad_Spill_Filter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Bad_Spill_Filter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Bad_Spill_Filter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(Bad_Spill_Filter);
