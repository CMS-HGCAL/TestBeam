//author: Thorben Quast
//data: 17th May
//structure representing reconstructed DATURA Telescope data. Reconstruction is performed using NtletTracking within the corryvreckan framework

#ifndef HGCalTBDATURATelescopeData_H
#define HGCalTBDATURATelescopeData_H

class HGCalTBDATURATelescopeData {
  public:  
    explicit HGCalTBDATURATelescopeData(int _track_ID): track_ID(_track_ID) {};
    HGCalTBDATURATelescopeData():track_ID(1) {};
    int track_ID;   

    void addPointForTracking(float _x, float _y, float _z, float _res_x, float _res_y) {
      x.push_back(_x);
      y.push_back(_y);
      z.push_back(_z);
      res_x.push_back(_res_x);
      res_y.push_back(_res_y);
      
      N_hitMultiplicity++;
    }
 

    //triplet tracking --> impact positions + associated chi2s of tracks...set in TelescopeProducer
    std::map<int, std::pair<double, double> > layerReferences;
    std::map<int, std::pair<double, double> > layerReferencesChi2;

    void addLayerReference(int l, double x, double y, double chi2_x, double chi2_y) {
        layerReferences[l] = std::make_pair(x, y);
        layerReferencesChi2[l] = std::make_pair(chi2_x, chi2_y);
    }

    std::pair<double, double> Extrapolation_XY(int l) const{
        if (layerReferences.find(l) != layerReferences.end())
            return layerReferences.at(l);
        else
            return std::make_pair(-999., -999.); 
    }  

    std::pair<double, double> Extrapolation_XY_Chi2(int l) const{
        if (layerReferencesChi2.find(l) != layerReferencesChi2.end())
            return layerReferencesChi2.at(l);
        else
            return std::make_pair(-999., -999.); 
    }   

  private:  
    int N_hitMultiplicity; 
    std::vector<float> x;    //mm
    std::vector<float> y;    //mm
    std::vector<float> z;    //mm, individual z values for each point due to alignment
    std::vector<float> res_x;    //resolution in mm
    std::vector<float> res_y;

    
};



#endif