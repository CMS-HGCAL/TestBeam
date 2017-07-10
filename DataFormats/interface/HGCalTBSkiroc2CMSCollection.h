#ifndef DATAFORMATS_HGCALTBSKIROC2CMSCOLLECTION_H
#define DATAFORMATS_HGCALTBSKIROC2CMSCOLLECTION_H

#include "DataFormats/Common/interface/EDCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBSkiroc2CMS.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

typedef edm::EDCollection<HGCalTBSkiroc2CMS> HGCalTBSkiroc2CMSCollection;
typedef edm::Ref<HGCalTBSkiroc2CMSCollection> HGCalTBSkiroc2CMSRef;
typedef edm::RefVector<HGCalTBSkiroc2CMSCollection> HGCalTBSkiroc2CMSRefs;
typedef edm::RefProd<HGCalTBSkiroc2CMSCollection> HGCalTBSkiroc2CMSRefProd;


#endif
