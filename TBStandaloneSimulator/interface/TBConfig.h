#ifndef TBCONFIG_H
#define TBCONFIG_H

#include <iostream>
#include <string>

struct TBConfig {
	TBConfig() {}
	TBConfig(std::string ctrlfile);
	~TBConfig() {}

	std::string macro;
	std::string geometry;
	bool savetracks;
};

std::ostream& operator<<(std::ostream& os, TBConfig& config);

#endif
