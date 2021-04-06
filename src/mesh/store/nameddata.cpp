
#include "nameddata.h"
#include "basis/utilities/packing.h"

using namespace espreso;

std::vector<std::string> NamedData::coordinateSuffixes = { "_X", "_Y", "_Z", "_XY", "_YZ", "_XZ", "_YX", "_ZY", "ZX" };
std::vector<std::string> NamedData::numberSuffixes = { "_1ST", "_2ND", "_3RD" };

NamedData::NamedData(int dimension, DataType datatype, const std::string &name)
: dimension(dimension), dataType(datatype), name(name), restriction(step::TYPE::TIME)
{

}

NamedData::NamedData(const char* &packedData)
{
	utils::unpack(dimension, packedData);
	utils::unpack(dataType, packedData);
	utils::unpack(name, packedData);
	utils::unpack(restriction, packedData);
	utils::unpack(data, packedData);
}

bool NamedData::withSuffixes() const
{
	return 1 < dimension;
}

bool NamedData::onlySuffixed() const
{
	return  dataType == DataType::TENSOR_ASYM ||
			dataType == DataType::TENSOR_SYMM ||
			(dataType == DataType::SCALAR && 1 < dimension);
}

int NamedData::nstatistics() const
{
	if (onlySuffixed() || dimension == 1) {
		return dimension;
	} else {
		return dimension + 1;
	}
}

std::string NamedData::suffix(int index) const
{
	switch (dataType) {
	case DataType::NUMBERED:
		return numberSuffixes[index];
	default:
		return coordinateSuffixes[index];
	}
}

size_t NamedData::packedSize()
{
	size_t packetSize = 0;
	packetSize += utils::packedSize(dimension);
	packetSize += utils::packedSize(dataType);
	packetSize += utils::packedSize(name);
	packetSize += utils::packedSize(restriction);
	packetSize += utils::packedSize(data);
	return packetSize;
}

void NamedData::pack(char *&p)
{
	utils::pack(dimension, p);
	utils::pack(dataType, p);
	utils::pack(name, p);
	utils::pack(restriction, p);
	utils::pack(data, p);
}
