#ifndef HGCalTBRunData_H
#define HGCalTBRunData_H
#include <string>
#include <map>

enum RUNTYPES {
  HGCAL_TB_PEDESTAL = 0,
  HGCAL_TB_BEAM , 
  HGCAL_TB_SIM,
  HGCAL_TB_SIM_FTF_BIC,
  HGCAL_TB_SIM_FTFP_BERT, 
  HGCAL_TB_SIM_FTFP_BERT_EML, 
  HGCAL_TB_SIM_FTFP_BERT_EMM,
  HGCAL_TB_SIM_FTFP_BERT_HP_EML,
  HGCAL_TB_SIM_FTFP_BERT_EMY,
  HGCAL_TB_SIM_QGSP_BERT,
  HGCAL_TB_SIM_QGSP_BERT_EML,
  HGCAL_TB_SIM_QGSP_BERT_HP_EML,
  HGCAL_TB_SIM_QGSP_FTFP_BERT,
  HGCAL_TB_SIM_QGSP_FTFP_BERT_EML,
  HGCAL_TB_SIM_QGSP_FTFP_BERT_EML_New,
  HGCAL_TB_SIM_QGSP_FTFP_BERT_EMM,
  HGCAL_TB_SIM_GAN,
  HGCAL_TB_SIM_GAN_1,
  HGCAL_TB_SIM_GAN_2,
  HGCAL_TB_SIM_GAN_3
};


template <class T>
class UserRecords { 
   private: 
      std::map<std::string, T> elements;

   public: 
      UserRecords<T>() {};

      void add(const std::string key, T const& elem) {
        elements[key] = elem;
      };   
      bool has(const std::string key) const{
        typename std::map<std::string, T>::const_iterator it = elements.find(key);
        return (it!=elements.end());
      };                
      T get(const std::string key) const{
        return this->elements.at(key);  //"[]" operator is a non-const labelled operator
      };              
}; 


class RunData {
  public:
    explicit RunData(int config, int r, double e, std::string rt):  run(r), configuration(config), energy(e) {};
    RunData(): run(0), configuration(0), energy(0), pdgID(0), runType(HGCAL_TB_BEAM){};
    int run;
    int trigger;
    int event;
    int configuration;
    double energy;
    int pdgID;
    RUNTYPES runType;
    
    UserRecords<bool> booleanUserRecords;
    UserRecords<double> doubleUserRecords;

};

typedef std::map<int, RunData> runMap;


#endif