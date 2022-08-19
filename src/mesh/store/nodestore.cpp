
#include "store.h"
#include "nodestore.h"
#include "statisticsstore.h"

#include "mesh/mesh.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

using namespace espreso;

NodeStore::NodeStore()
: size(0),
  distribution({0, 0}),
  IDs(NULL),
  elements(NULL),
  inputOffset(NULL),
  outputOffset(NULL),
  originCoordinates(NULL),
  coordinates(NULL),
  ranks(NULL),
  domains(NULL),
  eregions(NULL),
  bregions(NULL)
{

}

size_t NodeStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(size);
	packedSize += utils::packedSize(uniqInfo.nhalo);
	packedSize += utils::packedSize(uniqInfo.offset);
	packedSize += utils::packedSize(uniqInfo.size);
	packedSize += utils::packedSize(uniqInfo.totalSize);
	packedSize += utils::packedSize(uniqInfo.position);
	packedSize += utils::packedSize(distribution);

	packedSize += utils::packedSize(IDs);
	packedSize += utils::packedSize(elements);
	packedSize += utils::packedSize(inputOffset);
	packedSize += utils::packedSize(outputOffset);
	packedSize += utils::packedSize(originCoordinates);
	packedSize += utils::packedSize(coordinates);
	packedSize += utils::packedSize(ranks);
	packedSize += utils::packedSize(domains);
	packedSize += utils::packedSize(eregions);
	packedSize += utils::packedSize(bregions);

	packedSize += utils::packedSize(data.size());
	for (size_t i = 0; i < data.size(); i++) {
		packedSize += data[i]->packedSize();
	}

	return packedSize;
}

void NodeStore::packFull(char* &p) const
{
	utils::pack(size, p);
	utils::pack(uniqInfo.nhalo, p);
	utils::pack(uniqInfo.offset, p);
	utils::pack(uniqInfo.size, p);
	utils::pack(uniqInfo.totalSize, p);
	utils::pack(uniqInfo.position, p);
	utils::pack(distribution, p);

	utils::pack(IDs, p);
	utils::pack(elements, p);
	utils::pack(inputOffset, p);
	utils::pack(outputOffset, p);
	utils::pack(originCoordinates, p);
	utils::pack(coordinates, p);
	utils::pack(ranks, p);
	utils::pack(domains, p);
	utils::pack(eregions, p);
	utils::pack(bregions, p);

	utils::pack(data.size(), p);
	for (size_t i = 0; i < data.size(); i++) {
		data[i]->pack(p);
	}
}

void NodeStore::unpackFull(const char* &p)
{
	utils::unpack(size, p);
	utils::unpack(uniqInfo.nhalo, p);
	utils::unpack(uniqInfo.offset, p);
	utils::unpack(uniqInfo.size, p);
	utils::unpack(uniqInfo.totalSize, p);
	utils::unpack(uniqInfo.position, p);
	utils::unpack(distribution, p);

	utils::unpack(IDs, p);
	utils::unpack(elements, p);
	utils::unpack(inputOffset, p);
	utils::unpack(outputOffset, p);
	utils::unpack(originCoordinates, p);
	utils::unpack(coordinates, p);
	utils::unpack(ranks, p);
	utils::unpack(domains, p);
	utils::unpack(eregions, p);
	utils::unpack(bregions, p);

	size_t size;
	utils::unpack(size, p);
	for (size_t i = 0; i < size; i++) {
		data.push_back(new NodeData(p));
	}
}

size_t NodeStore::packedSize() const
{
	return
			utils::packedSize(size) +
			utils::packedSize(uniqInfo.nhalo) +
			utils::packedSize(uniqInfo.offset) +
			utils::packedSize(uniqInfo.size) +
			utils::packedSize(uniqInfo.totalSize) +
			utils::packedSize(uniqInfo.position) +
			IDs->packedSize() +
			coordinates->packedSize();
}

void NodeStore::pack(char* &p) const
{
	utils::pack(size, p);
	utils::pack(uniqInfo.nhalo, p);
	utils::pack(uniqInfo.offset, p);
	utils::pack(uniqInfo.size, p);
	utils::pack(uniqInfo.totalSize, p);
	utils::pack(uniqInfo.position, p);
	IDs->pack(p);
	coordinates->pack(p);
}

void NodeStore::unpack(const char* &p)
{
	if (IDs == NULL) {
		IDs = new serializededata<esint, esint>(1, tarray<esint>(1, 0));
	}
	if (coordinates == NULL) {
		coordinates = new serializededata<esint, Point>(1, tarray<Point>(1, 0));
	}

	utils::unpack(size, p);
	utils::unpack(uniqInfo.nhalo, p);
	utils::unpack(uniqInfo.offset, p);
	utils::unpack(uniqInfo.size, p);
	utils::unpack(uniqInfo.totalSize, p);
	utils::unpack(uniqInfo.position, p);
	IDs->unpack(p);
	coordinates->unpack(p);
}

size_t NodeStore::packedDataHeaderSize() const
{
	size_t size = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->name.size()) {
			size += utils::packedSize(data[i]->dimension);
			size += utils::packedSize(data[i]->dataType);
			size += utils::packedSize(data[i]->name);
			size += utils::packedSize(data[i]->restriction);
		}
	}
	return size;
}

void NodeStore::packDataHeader(char* &p) const
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->name.size()) {
			size += 1;
		}
	}
	utils::pack(size, p);
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->name.size()) {
			utils::pack(data[i]->dimension, p);
			utils::pack(data[i]->dataType, p);
			utils::pack(data[i]->name, p);
			utils::pack(data[i]->restriction, p);
		}
	}
}

void NodeStore::unpackDataHeader(const char* &p)
{
	size_t size;
	utils::unpack(size, p);
	for (size_t i = 0; i < size; i++) {
		data.push_back(new NodeData(0, NamedData::DataType::VECTOR, {}));
		utils::unpack(data[i]->dimension, p);
		utils::unpack(data[i]->dataType, p);
		utils::unpack(data[i]->name, p);
		utils::unpack(data[i]->restriction, p);
	}
}

size_t NodeStore::packedDataSize() const
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->name.size()) {
			size += utils::packedSize(data[i]->data);
		}
	}
	return size;
}

void NodeStore::packData(char* &p) const
{
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->name.size()) {
			utils::pack(data[i]->data, p);
		}
	}
}

void NodeStore::unpackData(const char* &p)
{
	for (size_t i = 0; i < data.size(); i++) {
		utils::unpack(data[i]->data, p);
	}
}

NodeStore::~NodeStore()
{
	if (IDs != NULL) { delete IDs; }
	if (elements != NULL) { delete elements; }

	if (inputOffset != NULL) { delete inputOffset; }
	if (outputOffset != NULL) { delete outputOffset; }
	if (originCoordinates != NULL) { delete originCoordinates; }
	if (coordinates != NULL) { delete coordinates; }
	if (ranks != NULL) { delete ranks; }
	if (domains != NULL) { delete domains; }
	if (eregions != NULL) { delete eregions; }
	if (bregions != NULL) { delete bregions; }

	for (size_t i = 0; i < data.size(); i++) {
		delete data[i];
	}
}

void NodeStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(info::mpi::rank) + ".txt");

	Store::storedata(os, "IDs", IDs);
	Store::storedata(os, "elements", elements);

	Store::storedata(os, "coordinates", coordinates);
	Store::storedata(os, "ranks", ranks);
	Store::storedata(os, "domains", domains);
	Store::storedata(os, "eregions", eregions);
	Store::storedata(os, "bregions", bregions);
}

void NodeStore::permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (IDs != NULL) { IDs->permute(permutation, distribution); }
	if (elements != NULL) { elements->permute(permutation, distribution); }

	if (inputOffset != NULL) { inputOffset->permute(permutation, distribution); }
	if (outputOffset != NULL) { outputOffset->permute(permutation, distribution); }
	if (originCoordinates != NULL) { originCoordinates->permute(permutation, distribution); }
	if (coordinates != NULL) { coordinates->permute(permutation, distribution); }
	if (ranks != NULL) { ranks->permute(permutation, distribution); }
	if (domains != NULL) { domains->permute(permutation, distribution); }
	if (eregions != NULL) { eregions->permute(permutation, distribution); }
	if (bregions != NULL) { bregions->permute(permutation, distribution); }
}

std::vector<esint> NodeStore::gatherNodeDistribution()
{
	return Store::gatherDistribution(size);
}

std::vector<esint> NodeStore::gatherUniqueNodeDistribution()
{
	return Store::gatherDistribution(uniqInfo.size);
}

NodeData* NodeStore::appendData(int dimension, NamedData::DataType datatype, const std::string &name, step::TYPE restriction)
{
	data.push_back(new NodeData(dimension, datatype, name));
	data.back()->restriction = restriction;
	data.back()->data.resize(dimension * size);
	return data.back();
}

void NodeData::statistics(const tarray<esint> &nodes, esint totalsize, Statistics *statistics) const
{
	for (int d = 0; d < nstatistics(); d++) {
		(statistics + d)->reset();
	}

	auto nranks = info::mesh->nodes->ranks->begin();
	esint prev = 0;
	esint doffset = nstatistics() != dimension ? 1 : 0;
	for (auto n = nodes.begin(); n != nodes.end(); prev = *n++) {
		nranks += *n - prev;
		if (*nranks->begin() == info::mpi::rank) {
			double value = 0;
			for (int d = 0; d < dimension; d++) {
				value += store[*n * dimension + d] * store[*n * dimension + d];
				(statistics + d + doffset)->min    = std::min((statistics + d + doffset)->min, store[*n * dimension + d]);
				(statistics + d + doffset)->max    = std::max((statistics + d + doffset)->max, store[*n * dimension + d]);
				(statistics + d + doffset)->avg   += store[*n * dimension + d];
				(statistics + d + doffset)->norm  += store[*n * dimension + d] * store[*n * dimension + d];
				(statistics + d + doffset)->absmin = std::min((statistics + d + doffset)->absmin, std::fabs(store[*n * dimension + d]));
				(statistics + d + doffset)->absmax = std::max((statistics + d + doffset)->absmax, std::fabs(store[*n * dimension + d]));
			}
			if (dataType == DataType::VECTOR) {
				value = std::sqrt(value);
				statistics->min = std::min(statistics->min, value);
				statistics->max = std::max(statistics->max, value);
				statistics->avg += value;
				statistics->norm += value * value;
				statistics->absmin = std::min(statistics->absmin, std::fabs(value));
				statistics->absmax = std::max(statistics->absmax, std::fabs(value));
			}
		}
	}

	std::vector<Statistics> global(nstatistics());
	Communication::allReduce(statistics, global.data(), nstatistics(), MPITools::operations->STATISTICS, MPITools::operations->mergeStatistics);
	memcpy(statistics, global.data(), sizeof(Statistics) * nstatistics());

	for (int i = 0; i < nstatistics(); i++) {
		(statistics + i)->avg /= totalsize;
		(statistics + i)->norm = std::sqrt((statistics + i)->norm);
	}
}


