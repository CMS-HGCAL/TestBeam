// -------------------------------------------------------------------------
// Configuration for standalone simulator
// -------------------------------------------------------------------------
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "HGCal/TBStandaloneSimulator/interface/TBConfig.h"
// -------------------------------------------------------------------------
using namespace std;

TBConfig::TBConfig(string ctrlfile)
	: macro(""),
	  geometry(""),
	  savetracks(false)
{
	ifstream fin(ctrlfile.c_str());
	if ( !fin.good() ) {
		cout << "** TConfig ** cannot open control file "
		     << ctrlfile << endl;
		exit(0);
	}

	string token, line;
	while (getline(fin, line)) {
		if ( line == "" ) continue;
		if ( line.substr(0, 1) == "#" ) continue;
		if ( line.substr(0, 1) == "%" ) continue;
		if ( line.substr(0, 2) == "//" )continue;

		istringstream sin(line);
		sin >> token >> line;
		if      ( token.substr(0, 3) == "mac" )
			macro = line;
		else if ( token.substr(0, 3) == "geo" )
			geometry = line;
		else if ( token.substr(0, 5) == "savet" ) {
			istringstream sinn(line);
			sinn >> savetracks;
		}
	}
	fin.close();
}

std::ostream& operator<<(std::ostream& os, TBConfig& config)
{
	os << "TBStandaloneSimulator configuration" << endl;
	os << "   macro      " << config.macro << endl;
	os << "   geometry   " << config.geometry << endl;
	os << "   savetracks " << config.savetracks;
	return os;
}
