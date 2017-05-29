
#include "solution.h"

#include "../mesh/settings/property.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/elementtypes.h"
#include "../basis/logging/logging.h"

using namespace espreso;

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, const std::vector<Property> &properties, const std::vector<std::vector<double> > &data)
: name(name), eType(eType), properties(properties), data(data), _offset(static_cast<int>(Property::SIZE), -1), _statistic(eType, mesh, data, properties.size())
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
}

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, const std::vector<Property> &properties)
: name(name), eType(eType), properties(properties), data(_data), _offset(static_cast<int>(Property::SIZE)), _statistic(eType, mesh, data, properties.size())
{
	for (size_t p = 0; p < properties.size(); p++) {
		_offset[static_cast<int>(properties[p])] = p;
	}
	_data.resize(mesh.parts());
	for (size_t p = 0; p < mesh.parts(); p++) {
		switch (eType) {
		case ElementType::ELEMENTS:
			_data[p].resize(properties.size() * (mesh.getPartition()[p + 1] - mesh.getPartition()[p]));
			break;
		case ElementType::NODES:
			_data[p].resize(properties.size() * mesh.coordinates().localSize(p));
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid solution element type.";
		}
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


