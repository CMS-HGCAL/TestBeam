// ------------------------------------------------------------------------
// File: TBGeometry
// Description: read test beam geometry from Python config file and
//              make available in a convenient format
// Created: 28-Apr-2016 Harrison B. Prosper
// ------------------------------------------------------------------------
#include <Python.h>
#include <iostream>
#include <cassert>
#include "HGCal/TBStandaloneSimulator/interface/TBGeometry.h"
#include "TPython.h"

using namespace std;

namespace
{
void setString(TBGeometry::Element& c, PyObject* o, string name)
{
	//cout << "\tsetString(" << name << ")" << endl;
	assert(o);
	PyObject* q = PyDict_GetItemString(o, name.c_str());
	assert(q);
	c.smap[name] = string(PyString_AsString(q));
}
void setDouble(TBGeometry::Element& c, PyObject* o, string name)
{
	//cout << "\tsetDouble(" << name << ")" << endl;
	assert(o);
	PyObject* q = PyDict_GetItemString(o, name.c_str());
	assert(q);
	c.dmap[name] = PyFloat_AsDouble(q);
}
void setInteger(TBGeometry::Element& c, PyObject* o, string name)
{
	//cout << "\tsetInteger(" << name << ")" << endl;
	assert(o);
	PyObject* q = PyDict_GetItemString(o, name.c_str());
	assert(q);
	c.imap[name] = static_cast<int>(PyInt_AsLong(q));
}

template <typename T>
void debug(map<string, T>& o)
{
	for(typename map<string, T>::iterator it = o.begin(); it != o.end(); it++)
		cout  << "  " << it->first << ":\t" << it->second << endl;
}
};


TBGeometry::TBGeometry(string modulename)
	: _modulename(modulename),
	  _world(TBGeometry::Element()),
	  _header(TBGeometry::Element()),
	  _geometry(std::vector<TBGeometry::Element>()),
	  _sensitive(std::map<int, int>())
{
	cout << endl
	     << "== TBGeometry =="
	     << endl << endl;

	// import Python function createGeometry
	char cmd[2048];
	sprintf(cmd,
	        "from HGCal.TBStandaloneSimulator.TBGeometryUtil import "
	        "createGeometry");
	cout << cmd << endl;
	TPython::Exec(cmd);

	// execute createGeometry
	sprintf(cmd, "createGeometry('%s')", _modulename.c_str());
	cout << endl << cmd << endl;
	PyObject* geom = (PyObject*)TPython::Eval(cmd);
	assert(geom);

	// -----------------------------------------------------
	// get world information
	// -----------------------------------------------------
	PyObject* world = PyDict_GetItemString(geom, "world");
	assert(world);
	setString(_world, world, "shape");
	setString(_world, world, "units");
	setDouble(_world, world, "xside");
	setDouble(_world, world, "yside");
	setDouble(_world, world, "zside");

	//debug(_world.smap);
	//debug(_world.dmap);


	// -----------------------------------------------------
	// get geometry header information
	// -----------------------------------------------------
	PyObject* header = PyDict_GetItemString(geom, "header");
	if ( header ) {
		setInteger(_header, header, "model");
		setInteger(_header, header, "version");
	} else {
		_header.imap["model"] = 1;
		_header.imap["version"] = 1;
	}
	// -----------------------------------------------------
	// get geometry element information
	// -----------------------------------------------------
	PyObject* geometry = PyDict_GetItemString(geom, "geometry");
	assert(geometry);

	// loop over elements (layers)
	int nelements = PyList_Size(geometry);
	for(int c = 0; c < nelements; c++) {
		PyObject* item = PyList_GetItem(geometry, c);
		PyObject* keys = PyDict_Keys(item);
		int nkeys = PyList_Size(keys);

		TBGeometry::Element element;
		setInteger(element, item, "sensitive");
		setString(element, item, "shape");
		setString(element, item, "material");
		setString(element, item, "units");
		for(int ii = 0; ii < nkeys; ii++) {
			PyObject* key = PyList_GetItem(keys, ii);
			string name = string(PyString_AsString(key));
			if ( name == "sensitive" ) continue;
			if ( name == "shape" ) continue;
			if ( name == "material" ) continue;
			if ( name == "units" ) continue;
			if ( element.dmap.find(name) != element.dmap.end() ) continue;
			if ( name == "first" || name == "last" )
				setInteger(element, item, name);
			else
				setDouble(element, item, name);
		}
		_geometry.push_back(element);

		//cout << endl;
		//debug(_geometry.back().smap);
		//debug(_geometry.back().dmap);
	}

	// -----------------------------------------------------
	// get indices of sensitive layers
	// -----------------------------------------------------
	int layer = 0;
	for(size_t c = 0; c < _geometry.size(); c++) {
		TBGeometry::Element& element = _geometry[c];
		if ( element.imap["sensitive"] ) {
			layer++;
			_sensitive[layer] = c;
		}
	}
}

TBGeometry::~TBGeometry()
{
}

//
TBGeometry::Element
TBGeometry::operator()(string item, int index)
{
	if      ( item == "world" ) {
		return _world;
	} else if ( item == "sensitive" ) {
		if ( _sensitive.find(index) != _sensitive.end() )
			return _geometry[_sensitive[index]];
		else
			return TBGeometry::Element();
	} else if ( item == "geometry" ) {
		if ( index < 0 )
			return TBGeometry::Element();
		if ( index >= (int)_geometry.size() )
			return TBGeometry::Element();

		return _geometry[index];
	} else
		return TBGeometry::Element();
}

std::ostream& operator<<(std::ostream& os, TBGeometry& o)
{
	char record[1024];
	os << "World" << endl << endl;
	TBGeometry::Element world = o("world");
	for(map<string, string>::iterator it = world.smap.begin();
	        it != world.smap.end();
	        it++) {
		sprintf(record, "  %-16s:\t%s",
		        (it->first).c_str(), (it->second).c_str());
		os << record << endl;
	}
	for(map<string, double>::iterator it = world.dmap.begin();
	        it != world.dmap.end();
	        it++) {
		sprintf(record, "  %-16s:\t%f",
		        (it->first).c_str(), it->second);
		os << record << endl;
	}

	os << endl;
	os << "Geometry - number of elements: " << o.size() << endl;

	for(size_t c = 0; c < o.size(); c++) {
		TBGeometry::Element element = o("geometry", c);
		os << endl;
		if ( element.imap["first"] )
			os << "BEGIN(SamplingSection)"
			   << endl;

		os << " element " << c;
		if ( element.imap["sensitive"] )
			os << "\t== sensitive ==";
		os << endl;

		for(map<string, string>::iterator it = element.smap.begin();
		        it != element.smap.end();
		        it++) {
			sprintf(record, "  %-16s:\t%s",
			        (it->first).c_str(), (it->second).c_str());
			os << record << endl;
		}

		for(map<string, double>::iterator it = element.dmap.begin();
		        it != element.dmap.end();
		        it++) {
			sprintf(record, "  %-16s:\t%f",
			        (it->first).c_str(), it->second);
			os << record << endl;
		}

		if ( element.imap["last"] )
			os << "END(SamplingSection)" << endl;
	}
	return os;
}
