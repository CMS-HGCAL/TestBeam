#ifndef HGCALTBSORTINGHELPER_HH
#define HGCALTBSORTINGHELPER_HH

template <typename T, typename S>
  class SortByLayer
{
 public : 
  SortByLayer(){;}
  ~SortByLayer(){;}
    
  static bool sort(T t, S s){ return t.layer() < s.layer(); }
};

template <typename T, typename S>
  class SortBySize
{
 public : 
  SortBySize(){;}
  ~SortBySize(){;}
    
  static bool sort(T t, S s){ return t.size() > s.size(); }
};

template <typename T, typename S>
  class SortByEnergy
{
 public : 
  SortByEnergy(){;}
  ~SortByEnergy(){;}
    
  static bool sort(T t, S s){ return t.energyHigh() > s.energyHigh(); }
};

#endif
