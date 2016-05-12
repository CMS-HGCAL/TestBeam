#ifndef TBGEOMETRY_H
#define TBGEOMETRY_H

#include <string>
#include <vector>
#include <map>

class TBGeometry
{
public:
	struct Element {
		std::map<std::string, std::string> smap;
		std::map<std::string, double> dmap;
		std::map<std::string, int> imap;
	};

	TBGeometry() {}
	TBGeometry(std::string modulename);
	virtual ~TBGeometry();

	size_t size()
	{
		return _geometry.size();
	}
	int    model()
	{
		return _header.imap["model"];
	}
	int    version()
	{
		return _header.imap["version"];
	}

	std::string material(std::string key)
	{
		if ( _header.smap.find(key) != _header.smap.end() )
			return _header.smap[key];
		else
			return "";
	}

	TBGeometry::Element operator()(std::string, int index = 0);


private:
	std::string                      _modulename;
	TBGeometry::Element              _world;
	TBGeometry::Element              _header;
	std::vector<TBGeometry::Element> _geometry;
	std::map<int, int>               _sensitive;
};

std::ostream& operator<<(std::ostream&, TBGeometry&);

#endif
