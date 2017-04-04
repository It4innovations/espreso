
#include "solution.h"

#include "../mesh/settings/property.h"

using namespace espreso;

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, const std::vector<Property> &properties, const std::vector<std::vector<double> > &data)
: name(name), eType(eType), properties(properties.size()), data(data), _offset(static_cast<int>(Property::SIZE), -1), _statistic(eType, mesh, data)
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
}

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, const std::vector<Property> &properties)
: name(name), eType(eType), properties(properties.size()), data(_data), _offset(static_cast<int>(Property::SIZE)), _statistic(eType, mesh, data)
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
}

bool Solution::hasProperty(Property property) const
{
	return _offset[static_cast<int>(property)] != -1;
}

void Solution::computeStatisticalData()
{
	_statistic.compute();
}

double Solution::getStatisticalData(Property property, StatisticalData data, const Region *region) const
{
	return _statistic.get(region, _offset[static_cast<int>(property)], data);
}


