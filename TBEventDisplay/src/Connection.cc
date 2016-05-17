///////////////////////////////////////////////////////////////////////
// File:    Connection.cc
// Purpose: Model a connection between a Signal and a Slot
// Created: Summer-2002 Harrison B. Prosper
// Updated: 05-Jun-2008 HBP Adapt to CMS
///////////////////////////////////////////////////////////////////////
//$Revision: 1.2 $

//#include "PhysicsTools/PyGui/interface/Connection.h"
#include "HGCal/TBEventDisplay/interface/Connection.h"
#include <iostream>
#include <iomanip>

using namespace std;

Connection::Connection(TQObject *sender,   const char *signal,
                       PyObject *receiver, const char *method)
	: _sender  (sender),
	  _signal  (signal)
{
	_slot = new Slot(receiver, method);

	if ( noArgs(_signal) )
		TQObject::Connect(_sender, _signal.c_str(),
		                  "Slot", _slot, "handleSignal()");

	else
		TQObject::Connect(_sender, _signal.c_str(),
		                  "Slot", _slot, "handleSignal(int)");
}

Connection::~Connection()
{
	if ( noArgs(_signal) )
		TQObject::Disconnect(_sender, _signal.c_str(),
		                     _slot, "handleSignal()");
	else
		TQObject::Disconnect(_sender, _signal.c_str(),
		                     _slot, "handleSignal(int)");
	delete _slot;
	_slot = 0;
}

Bool_t Connection::noArgs(string& signal)
{
	return signal.find("()") < signal.length();
}
