
#include "property.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/logging/logging.h"

namespace espreso {

std::istream& operator>>(std::istream& is, Property &property)
{
	std::string name;
	is >> name;
	for (size_t p = 0; p < (size_t)Property::SIZE; p++) {
		std::stringstream pname;
		pname << (Property)p;
		if (StringCompare::caseInsensitiveEq(name, pname.str())) {
			property = (Property)p;
			return is;
		}
	}

	ESINFO(ERROR) << "Unknown property '" << name << "'\n";
	return is;
}

}



