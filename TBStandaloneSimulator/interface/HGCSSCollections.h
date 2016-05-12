#ifndef HGCSSCOLLECTIONS_H
#define HGCSSCOLLECTIONS_H

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "HGCal/TBStandaloneSimulator/interface/HGCSSSimHit.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSGenParticle.hh"
#include "HGCal/TBStandaloneSimulator/interface/HGCSSSamplingSection.hh"

typedef edm::SortedCollection<HGCSSSimHit> HGCSSSimHitCollection;
typedef edm::Ref<HGCSSSimHitCollection> HGCSSSimHitRef;
typedef edm::RefVector<HGCSSSimHitCollection> HGCSSSimHitRefVector;
typedef edm::RefProd<HGCSSSimHitCollection> HGCSSSimHitRefProd;

typedef edm::SortedCollection<HGCSSGenParticle> HGCSSGenParticleCollection;
typedef edm::Ref<HGCSSGenParticleCollection> HGCSSGenParticleRef;
typedef edm::RefVector<HGCSSGenParticleCollection> HGCSSGenParticleRefVector;
typedef edm::RefProd<HGCSSGenParticleCollection> HGCSSGenParticleRefProd;

typedef edm::SortedCollection<HGCSSSamplingSection>
HGCSSSamplingSectionCollection;
typedef edm::Ref<HGCSSSamplingSectionCollection>
HGCSSSamplingSectionRef;
typedef edm::RefVector<HGCSSSamplingSectionCollection>
HGCSSSamplingSectionRefVector;
typedef edm::RefProd<HGCSSSamplingSectionCollection>
HGCSSSamplingSectionRefProd;

#endif
