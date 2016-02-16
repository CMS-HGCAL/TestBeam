#include <stdio.h>
#include <string.h>
#include <sstream>
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

static HGCalTBDetId tb_detid_load(const char* buffer, int& ptr) {
  int v0,v1,v2,v3,v4,v5;
  int found=sscanf(buffer,"%d %d %d %d %d %d%n ",&v0,&v1,&v2,&v3,&v4,&v5,&ptr);
  if (found==6) {
    return HGCalTBDetId(abs(v1),v2,v3,v4,v5!=0);
  } else return HGCalTBDetId(0);  
}
static void tb_detid_store(HGCalTBDetId id, FILE* f) {
  fprintf(f,"%08x %4d %4d %4d %4d %d ",id.rawId(),id.layer()*id.zside(),id.sensor(),id.ix(),id.iy(),(id.isCalib()?(1):(0)));
}

bool HGCalCondObjectTextIO::load(const std::string& filename, HGCalCondObjectContainer<double>& cont) {
  FILE* f=fopen(filename.c_str(),"r");
  if (f==0) {
    fprintf(stderr,"Unable to open '%s'\n",filename.c_str());
    return false;
  }

  char buffer[100];
  // first line is the scheme code
  uint64_t code(0);
  buffer[0]=0; fgets(buffer,100,f);
  if (sscanf(buffer,"SCHEME_CODE %lu",&code)!=1) {
    fprintf(stderr,"Expected 'SCHEME_CODE <value>' on first line of file\n");
    fclose(f);
    return false;
  }
  // clear the container
  cont=HGCalCondObjectContainer<double>(cont.getNumberingScheme(),code);

  while (!feof(f)) {
    buffer[0]=0; fgets(buffer,100,f);
    // trim comments
    char* p_comment=index(buffer,'#');
    if (p_comment!=0) *p_comment=0;
    // try to unpack...
    int ptr;
    HGCalTBDetId id=tb_detid_load(buffer,ptr);
    if (!id.null()) {
      double value=atof(buffer+ptr);
      cont.set(id,value);
    }    
  }
  fclose(f);
  return true;
}
bool HGCalCondObjectTextIO::store(const std::string& filename, const HGCalCondObjectContainer<double>& cont) {
  FILE* f=fopen(filename.c_str(),"r");
  if (f==0) return false;

  fprintf(f,"SCHEME_CODE %lu\n",cont.schemeCode());
  fprintf(f,"# HEX  LAYER SENSOR  IX  IY  CAL  VALUE\n");
  for (size_t i=0; i<cont.size(); i++) {
    if (cont.get(i).id.null()) continue;
    tb_detid_store(HGCalTBDetId(cont.get(i).id),f);
    fprintf(f,"%.6g\n",cont.get(i).value);
  }
  fclose(f);
  return true;
}
