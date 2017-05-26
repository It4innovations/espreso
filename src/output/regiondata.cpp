
#include <cstring>

#include "regiondata.h"

#include "../assembler/solution.h"
#include "../mesh/structures/elementtypes.h"

using namespace espreso::output;

DataArrays::~DataArrays()
{
	clear();
}

void DataArrays::clear()
{
	for (auto it = elementDataDouble.begin(); it != elementDataDouble.end(); ++it) {
		delete it->second.second;
	}
	for (auto it = elementDataInteger.begin(); it != elementDataInteger.end(); ++it) {
		delete it->second.second;
	}
	for (auto it = pointDataDouble.begin(); it != pointDataDouble.end(); ++it) {
		delete it->second.second;
	}
	for (auto it = pointDataInteger.begin(); it != pointDataInteger.end(); ++it) {
		delete it->second.second;
	}

	elementDataDouble.clear();
	elementDataInteger.clear();
	pointDataDouble.clear();
	pointDataInteger.clear();
}

std::vector<std::string> RegionData::pointDataNames() const
{
	std::vector<std::string> names;
	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		names.push_back(it->first);
	}
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		names.push_back(it->first);
	}
	return names;
}

std::vector<std::string> RegionData::cellDataNames() const
{
	std::vector<std::string> names;
	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		names.push_back(it->first);
	}
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		names.push_back(it->first);
	}
	return names;
}

void RegionData::clearData()
{
	data.clear();
}

void RegionData::pack(std::vector<char> &data) const
{

}

void RegionData::unpack(char* &data)
{

}



