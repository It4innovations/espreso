
#include "solution.h"

#include "../mesh/settings/property.h"

using namespace espreso;

Solution::Solution(const std::string &name, store::ElementType eType, const std::vector<Property> &properties, const std::vector<std::vector<double> > &data)
: name(name), eType(eType), properties(properties.size()), data(data), _offset(static_cast<int>(Property::SIZE))
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
}

Solution::Solution(const std::string &name, store::ElementType eType, const std::vector<Property> &properties)
: name(name), eType(eType), properties(properties.size()), data(_data), _offset(static_cast<int>(Property::SIZE))
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
}


