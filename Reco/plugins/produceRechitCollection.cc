#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include <iostream>
#include "TMath.h"
#include "TRandom3.h"

using namespace std;

class produceRechitCollection : public edm::EDProducer{

      public:
             produceRechitCollection(const edm::ParameterSet&);
             virtual void produce(edm::Event&, const edm::EventSetup&);             
      private:
              std::string outputCollName1;     ///<label name of collection made by this producer
              TRandom3 Rand1,Rand2;
    };

produceRechitCollection::produceRechitCollection(const edm::ParameterSet& cfg)
       :outputCollName1(cfg.getParameter<std::string>("OutputCollectionName1"))
{
   Rand1.SetSeed(15234525);
   Rand2.SetSeed(34252345);
   produces <HGCalTBRecHitCollection>(outputCollName1);
  }

void produceRechitCollection::produce(edm::Event& event, const edm::EventSetup&){

    auto_ptr<HGCalTBRecHitCollection> Rechits(new HGCalTBRecHitCollection);
    int Layer=1;
    int Sensor=1;
    unsigned int mode = 0;
    for(int i=-7;i<8;i++){
       for(int j=-7;j<8;j++){
           HGCalTBDetId Id(Layer,Sensor,i,j,false);
//           std::cout << std::endl<<HGCalTBRecHit(Id, 65 - fabs(i*j),mode) << std::endl;
           Rechits->push_back(HGCalTBRecHit(Id,  Rand2.Gaus(4,1) + fabs(TMath::Gaus(sqrt(i*i + j*j),0.,5.)*Rand1.Gaus(10.,5.)),0,mode));  
         }
      } 
    event.put(Rechits, outputCollName1);
  }
// Should there be a destructor ?? 
DEFINE_FWK_MODULE(produceRechitCollection);
