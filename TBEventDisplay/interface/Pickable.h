#ifndef PICKABLE_H
#define PICKABLE_H
//-----------------------------------------------------------------------------
// Package:    PhysicsTools/PyGui
// Created:    Created: April 2016 HBP
//-----------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <string>
//-----------------------------------------------------------------------------
#include "TQObject.h"
#include "TEveSelection.h"
//-----------------------------------------------------------------------------
class Pickable : public TQObject
{
public:
	Pickable();
	virtual ~Pickable();

	void AddElement(TEveElement* element);
	void Selected(TEveElement* element);
	void Clear();
	TEveElement* operator[](int id);

	void Selected(int id);     //*SIGNAL*
	void Cleared();            //*SIGNAL*

private:
	TEveSelection* _selection;
	std::map<TEveElement*, int> _element2id;
	std::map<int, TEveElement*> _id2element;
	int _id;

public:
	ClassDef(Pickable, 0)      // Needed to get signals to work
};

#endif
