
#include "meshinfo.h"

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

void MeshInfo::clearData()
{
	data.clear();
	solutions.clear();
}


