
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

template <typename TType>
static size_t packedVectorSize(const std::vector<TType> &vector)
{
	return sizeof(size_t) + sizeof(TType) * vector.size();
}

template <typename TType>
static size_t packedDataSize(const std::map<std::string, std::pair<size_t, std::vector<TType>* > > &data)
{
	size_t size = 0;
	size += sizeof(size_t);
	for (auto it = data.begin(); it != data.end(); ++it) {
		size += sizeof(size_t) + it->first.size();
		size += sizeof(size_t) + sizeof(size_t) + sizeof(TType) * it->second.second->size();
	}
	return size;
}

template <typename TType>
static void packVector(char* &p, std::vector<char> &packedData, const std::vector<TType> &vector)
{
	size_t datasize = vector.size();
	memcpy(p, &datasize, sizeof(size_t));
	p += sizeof(size_t);
	memcpy(p, vector.data(), sizeof(TType) * vector.size());
	p += sizeof(TType) * vector.size();
}

template <typename TType>
static void packData(char* &p, std::vector<char> &packedData, const std::map<std::string, std::pair<size_t, std::vector<TType>* > > &data)
{
	size_t datasize = data.size();
	memcpy(p, &datasize, sizeof(size_t));
	p += sizeof(size_t);

	for (auto it = data.begin(); it != data.end(); ++it) {
		datasize = it->first.size();
		memcpy(p, &datasize, sizeof(size_t));
		p += sizeof(size_t);
		memcpy(p, it->first.c_str(), it->first.size());
		p += it->first.size();

		datasize = it->second.first;
		memcpy(p, &datasize, sizeof(size_t));
		p += sizeof(size_t);

		datasize = it->second.second->size();
		memcpy(p, &datasize, sizeof(size_t));
		p += sizeof(size_t);
		memcpy(p, it->second.second->data(), sizeof(TType) * it->second.second->size());
		p += sizeof(TType) * it->second.second->size();
	}
}

template <typename TType>
static void unpackVector(char* &p, std::vector<TType> &vector)
{
	size_t datasize;
	memcpy(&datasize, p, sizeof(size_t));
	p += sizeof(size_t);
	vector.resize(datasize);
	memcpy(vector.data(), p, sizeof(TType) * vector.size());
	p += sizeof(TType) * vector.size();
}

template <typename TType>
static void unpackData(char* &p, std::map<std::string, std::pair<size_t, std::vector<TType>* > > &data)
{
	size_t vectors;
	memcpy(&vectors, p, sizeof(size_t));
	p += sizeof(size_t);

	for (size_t i = 0; i < vectors; i++) {
		size_t stringsize, properties, datasize;
		memcpy(&stringsize, p, sizeof(size_t));
		p += sizeof(size_t);
		std::string name(p, p + stringsize);
		p += stringsize;

		memcpy(&properties, p, sizeof(size_t));
		p += sizeof(size_t);

		memcpy(&datasize, p, sizeof(size_t));
		p += sizeof(size_t);

		std::vector<TType> *vector = new std::vector<TType>(datasize);
		memcpy(vector->data(), p, sizeof(TType) * vector->size());
		p += sizeof(TType) * vector->size();
		data[name] = std::make_pair(properties, vector);
	}
}

void RegionData::pack(std::vector<char> &data) const
{
	size_t size = 0;

	size += packedVectorSize(coordinates);
	size += packedVectorSize(elementsTypes);
	size += packedVectorSize(elementsNodes);
	size += packedVectorSize(elements);

	size += packedDataSize(this->data.pointDataDouble);
	size += packedDataSize(this->data.pointDataInteger);
	size += packedDataSize(this->data.elementDataDouble);
	size += packedDataSize(this->data.elementDataInteger);

	data.resize(data.size() + size);
	char *p = data.data() + data.size() - size;

	packVector(p, data, coordinates);
	packVector(p, data, elementsTypes);
	packVector(p, data, elementsNodes);
	packVector(p, data, elements);

	packData(p, data, this->data.pointDataDouble);
	packData(p, data, this->data.pointDataInteger);
	packData(p, data, this->data.elementDataDouble);
	packData(p, data, this->data.elementDataInteger);
}

void RegionData::unpack(char* &data)
{
	coordinates.clear();
	elementsTypes.clear();
	elementsNodes.clear();
	elements.clear();

	unpackVector(data, coordinates);
	unpackVector(data, elementsTypes);
	unpackVector(data, elementsNodes);
	unpackVector(data, elements);

	this->data.pointDataDouble.clear();
	this->data.pointDataInteger.clear();
	this->data.elementDataDouble.clear();
	this->data.elementDataInteger.clear();

	unpackData(data, this->data.pointDataDouble);
	unpackData(data, this->data.pointDataInteger);
	unpackData(data, this->data.elementDataDouble);
	unpackData(data, this->data.elementDataInteger);
}



