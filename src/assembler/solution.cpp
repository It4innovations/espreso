
#include "solution.h"

#include "../mesh/settings/property.h"

using namespace espreso;

Solution::Solution(const std::vector<Property> &properties, const std::vector<std::vector<double> > &data)
: data(data), properties(properties.size()), _offset(static_cast<int>(Property::SIZE))
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
}

Solution::Solution(const std::vector<Property> &properties)
: data(_data), properties(properties.size()), _offset(static_cast<int>(Property::SIZE))
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
}


