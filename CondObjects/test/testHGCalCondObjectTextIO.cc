#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/CondObjects/interface/HGCalTBNumberingScheme.h"
#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>

int main(int argc, char *argv[])
{
	std::string mapfile_in, valuefile_in;
	int opt;

	while ((opt = getopt(argc, argv, "m:v:")) != -1) {
		switch (opt) {
		case 'm':
			mapfile_in = optarg;
			break;
		case 'v':
			valuefile_in = optarg;
			break;
		default: /* '?' */
			fprintf(stderr, "Usage: %s [-m mapfile] [-v single_value_file]\n",
			        argv[0]);
			exit(EXIT_FAILURE);
		}
	}

	HGCalCondObjectTextIO io(HGCalTBNumberingScheme::scheme());

	if (!mapfile_in.empty()) {
		HGCalElectronicsMap emap;
		printf("Reading electronics map from '%s'...", mapfile_in.c_str());
		fflush(stdout);
		bool ok = io.load(mapfile_in, emap);
		if (ok) printf("OK\n");
		else printf("ERROR\n");
		printf("Rewriting electronics map to '%s.out'...", mapfile_in.c_str());
		fflush(stdout);
		ok = io.store(mapfile_in + ".out", emap);
		if (ok) printf("OK\n");
		else printf("ERROR\n");
	}

	if (!valuefile_in.empty()) {
		HGCalCondObjectContainer<float> values;
		printf("Reading single values from '%s'...", valuefile_in.c_str());
		fflush(stdout);
		bool ok = io.load(valuefile_in, values);
		if (ok) printf("OK\n");
		else printf("ERROR\n");
		printf("Rewriting single values to '%s.out'...", valuefile_in.c_str());
		fflush(stdout);
		ok = io.store(valuefile_in + ".out", values);
		if (ok) printf("OK\n");
		else printf("ERROR\n");
	}


	return 0;
}
