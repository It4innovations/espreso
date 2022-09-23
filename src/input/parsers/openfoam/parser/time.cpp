
#include "time.h"

#include "basis/utilities/parser.h"

#include <sstream>

using namespace espreso;

double OpenFOAMTime::value(const std::string &name)
{
	std::ifstream is(name);
	FoamFileHeader header; header.read(is);

	std::string parameter, value;
	while (!StringCompare::caseInsensitiveEq("value", parameter)) {
		is >> parameter >> value;
	}
	return std::stof(value);
}

