#ifndef SLOT_H
#define SLOT_H
/////////////////////////////////////////////////////////////////////////
// File:    Slot.h
// Purpose: Model a Slot for use in Signal/Slot communication.
// Created: Summer-2002 Harrison B. Prosper
// Updated: 05-Jun-2008 HBP Adapt to CMS
//          14-Apr-2011 HBP changed unsigned long
/////////////////////////////////////////////////////////////////////////
//$Revision: 1.2 $

#include "TQObject.h"
#include <string>
#include <vector>
#include <Python.h>

/**
 */
class Slot : public TQObject
{
private:
	PyObject*   _object;
	std::string _mstr;
	std::vector<char> _method;

public:

	/** RootCint requires a default constructor
	 */

	Slot();

	/**
	 */
	Slot(PyObject* object, const char* method);

	/**
	 */
	~Slot();

	/**
	*/
	void handleSignal(int id);

	/**
	 */
	void handleSignal();

	PyObject*  receiver() const
	{
		return _object;
	}

	const char*    method()   const
	{
		return _mstr.c_str();
	}

};

#endif
