//-----------------------------------------------------------------------------
/** PhysicsTools/PyGui/src/Pickable.cc

 Description: Wrapper around TEveSelection for use with Connection and Slot

 Implementation:
     As simple as possible
*/
// Created: April 2016 Harrison B. Prosper
//-----------------------------------------------------------------------------
#include <cassert>
#include <algorithm>
//#include "PhysicsTools/PyGui/interface/Pickable.h"
#include "HGCal/TBEventDisplay/interface/Pickable.h"
#include "TEveManager.h"
using namespace std;
//-----------------------------------------------------------------------------

ClassImp(Pickable); // Needed to get signals to work

Pickable::Pickable()
	: _selection(gEve->GetSelection()),
	  _element2id(map<TEveElement *, int>()),
	  _id2element(map<int, TEveElement * >()),
	  _id(0)
{
	_selection->Connect("SelectionAdded(TEveElement*)", "Pickable",
	                    this, "Selected(TEveElement*)");

	_selection->Connect("SelectionCleared()", "Pickable",
	                    this, "Cleared()");
}

Pickable::~Pickable()
{
	_selection->Disconnect("SelectionAdded(TEveElement*)",
	                       this, "Selected(TEveElement*)");

	_selection->Disconnect("SelectionCleared()",
	                       this, "Cleared()");
}

void
Pickable::AddElement(TEveElement* element)
{
	_element2id[element] = _id;
	_id2element[_id] = element;
	_id++;
}

void
Pickable::Clear()
{
	_selection->RemoveElements();
	_element2id.clear();
	_id2element.clear();
}

TEveElement*
Pickable::operator[](int id)
{
	TEveElement* element = 0;
	try {
		element = _id2element[id];
	} catch (...) {
		assert(element);
	}
	return element;
}

void
Pickable::Selected(int id)
{
	this->Emit("Selected(int)", id);
}

void
Pickable::Cleared()
{
	this->Emit("Cleared()");
}

void
Pickable::Selected(TEveElement* element)
{
	int id = -1;
	try {
		int id = _element2id[element];
		Selected(id);
	} catch (...) {
		assert(id > -1);
	}
}
