#include <stdio.h>
#include <string.h>
#include <sstream>
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

static HGCalTBDetId tb_detid_load(const char* buffer, int& ptr)
{
	int v1, v2, v3, v4, v5, v6;
	int found = sscanf(buffer, "%d %d %d %d %d %d %n ", &v1, &v2, &v3, &v4, &v5, &v6, &ptr);
	if (found == 6) {
		return HGCalTBDetId(abs(v1), v2, v3, v4, v5, v6);
	} else return HGCalTBDetId(0);
}
static void tb_detid_store(HGCalTBDetId id, FILE* f)
{
	fprintf(f, "%08x %5d %9d %9d %3d %3d %4d ", id.rawId(), id.layer()*id.zside(), id.sensorIU(), id.sensorIV(), id.iu(), id.iv(), id.cellType());
}

bool HGCalCondObjectTextIO::load(const std::string& filename, HGCalCondObjectContainer<float>& cont)
{
	FILE* f = fopen(filename.c_str(), "r");
	if (f == 0) {
		fprintf(stderr, "Unable to open '%s'\n", filename.c_str());
		return false;
	}

	char buffer[100];
	// first line is the scheme code
	uint64_t code(0);
	buffer[0] = 0;
	fgets(buffer, 100, f);
	if (sscanf(buffer, "SCHEME_CODE %lu", &code) != 1) {
		fprintf(stderr, "Expected 'SCHEME_CODE <value>' on first line of file\n");
		fclose(f);
		return false;
	}
	// clear the container
	cont = HGCalCondObjectContainer<float>(p_scheme, code);

	while (!feof(f)) {
		buffer[0] = 0;
		fgets(buffer, 100, f);
		// trim comments
		char* p_comment = index(buffer, '#');
		if (p_comment != 0) *p_comment = 0;
		// try to unpack...
		int ptr, dummy;
		// remove the first column
		const char* process = buffer;
		int found = sscanf(process, "%d %n ", &dummy, &ptr);
		if (found == 1) {
			process += ptr;
		} else continue;
		// get the detid
		HGCalTBDetId id = tb_detid_load(process, ptr);
		if (!id.null()) {
			process += ptr;
			float value = atof(process);
			cont.set(id, value);
		}
	}
	fclose(f);
	return true;
}
bool HGCalCondObjectTextIO::store(const std::string& filename, const HGCalCondObjectContainer<float>& cont)
{
	FILE* f = fopen(filename.c_str(), "w");
	if (f == 0) return false;

	fprintf(f, "SCHEME_CODE %lu\n", cont.schemeCode());
	fprintf(f, "#   CODE LAYER SENSOR_IU SENSOR_IV  IU  IV TYPE  VALUE\n");
	for (size_t i = 0; i < cont.size(); i++) {
		if (cont.get(i).id.null()) continue;
		tb_detid_store(HGCalTBDetId(cont.get(i).id), f);
		fprintf(f, " %.6g\n", cont.get(i).value);
	}
	fclose(f);
	return true;
}


bool HGCalCondObjectTextIO::load(const std::string& filename, HGCalElectronicsMap& emap)
{
	FILE* f = fopen(filename.c_str(), "r");
	if (f == 0) {
		fprintf(stderr, "Unable to open '%s'\n", filename.c_str());
		return false;
	}

	char buffer[100];
	while (!feof(f)) {
		buffer[0] = 0;
		fgets(buffer, 100, f);
		// trim comments
		char* p_comment = index(buffer, '#');
		if (p_comment != 0) *p_comment = 0;
		// try to unpack...
		int ptr;

		// electronics id sections
		int iskiroc, icell;
		const char* process = buffer;
		int found = sscanf(process, "%d %d %n", &iskiroc, &icell, &ptr);
		if (found == 2) {
			process += ptr;
		} else continue;
		// get the detid
		HGCalTBDetId id = tb_detid_load(process, ptr);
		if (!id.null()) {
			HGCalTBElectronicsId eid(iskiroc, icell);
			emap.insert(eid.rawId(), id);
		}
	}
	fclose(f);
	return true;
}
bool HGCalCondObjectTextIO::store(const std::string& filename, const HGCalElectronicsMap& emap)
{
	FILE* f = fopen(filename.c_str(), "w");
	if (f == 0) return false;

	fprintf(f, "# SKIROC CHAN | LAYER SENSOR_IU SENSOR_IV  IU  IV TYPE \n");
	for (size_t i = 0; i < emap.size(); i++) {
		HGCalTBElectronicsId eid(emap.eidAt(i));
		HGCalTBDetId did(emap.didAt(i));
		fprintf(f, "  %6d %4d   %5d %9d %9d %3d %3d %4d\n", eid.iskiroc(), eid.ichan(),
		        did.layer()*did.zside(), did.sensorIU(), did.sensorIV(), did.iu(), did.iv(), did.cellType());
	}


	fclose(f);
	return true;
}

bool HGCalCondObjectTextIO::load(const std::string& filename, HGCalTBDetectorLayout &hgcalGeom)
{
  FILE* f = fopen(filename.c_str(),"r");
  if (f == 0) {
    fprintf(stderr, "Unable to open '%s'\n", filename.c_str());
    return false;
  }
  int layerId,x,y,moduleId;
  float z;
  char det[10];
  char buffer[100];
  std::vector<HGCalTBLayer> layers;
  std::vector<HGCalTBLayer>::iterator it_layers;
  std::vector<HGCalTBModule> modules;
  std::vector<HGCalTBModule>::iterator it_modules;
  while(!feof(f)){
    buffer[0]=0;
    fgets(buffer,100,f);
    char* p_comment = index(buffer, '#');
    if (p_comment != 0) continue;
    int ptr=0;
    const char* process = buffer;
    int found = sscanf(process, "%d %f %s %d %d %d %n", &layerId, &z, det, &x, &y, &moduleId, &ptr);
    if (found == 6) {
      process += ptr;
      HGCalTBModule module(layerId, std::string(det), x, y, moduleId);
      HGCalTBLayer layer(layerId, z, std::string(det) );
      it_layers=std::find(layers.begin(), layers.end(), layer);
      modules.push_back(module);
      if( it_layers!=layers.end() )
	(*it_layers).add(module);
      else{
	layer.add(module);
	layers.push_back(layer);
      }
    } else continue;
  }
  for( it_layers=layers.begin(); it_layers!=layers.end(); ++it_layers ){
    hgcalGeom.add(*it_layers);
  }
  for( it_modules=modules.begin(); it_modules!=modules.end(); ++it_modules ){
    hgcalGeom.add(*it_modules);
  }
  fclose(f);
  return true;
}

bool HGCalCondObjectTextIO::load(const std::string& filename, HGCalTBADCConversionsMap &adcConvMap)
{
  FILE* f = fopen(filename.c_str(),"r");
  if (f == 0) {
    fprintf(stderr, "Unable to open '%s'\n", filename.c_str());
    return false;
  }
  int moduleId,asicId;
  float lg_hg_tr,lg_hg_cv,tot_lg_tr,tot_lg_cv,adc_to_mip;
  char buffer[100];
  while(!feof(f)){
    buffer[0]=0;
    fgets(buffer,100,f);
    char* p_comment = index(buffer, '#');
    if (p_comment != 0) continue;
    int ptr=0;
    const char* process = buffer;
    int found = sscanf(process, "%d %d %f %f %f %f %f %n", &moduleId, &asicId, &adc_to_mip,&lg_hg_tr, &lg_hg_cv,&tot_lg_tr, &tot_lg_cv, &ptr);
    if (found == 7) {
      process += ptr;
      ASIC_ADC_Conversions adcConv(moduleId,asicId,adc_to_mip,
				   lg_hg_tr,lg_hg_cv,
				   tot_lg_tr, tot_lg_cv);
      adcConvMap.addEntry( adcConv);
    } else continue;
  }
  fclose(f);
  return true;
}
