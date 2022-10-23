
#include "store.h"
#include "elementstore.h"
#include "statisticsstore.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"

#include "mesh/mesh.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

using namespace espreso;

ElementStore::ElementStore()
: offset(NULL),
  inputOffset(NULL),
  nodes(NULL),
  centers(NULL),

  body(NULL),
  contact(NULL),
  material(NULL),
  regions(NULL),
  epointers(NULL),

  faceNeighbors(NULL),
  edgeNeighbors(NULL),

  volumeIndices(NULL),

  stiffness(NULL)
{

}

size_t ElementStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(distribution.threads);
	packedSize += utils::packedSize(distribution.process);
	packedSize += utils::packedSize(distribution.code);

	packedSize += utils::packedSize(offset);
	packedSize += utils::packedSize(inputOffset);
	packedSize += utils::packedSize(nodes);
	packedSize += utils::packedSize(centers);
	packedSize += utils::packedSize(body);
	packedSize += utils::packedSize(contact);
	packedSize += utils::packedSize(material);
	packedSize += utils::packedSize(regions);
	packedSize += utils::packedSize(faceNeighbors);
	packedSize += utils::packedSize(edgeNeighbors);
	packedSize += utils::packedSize(volumeIndices);
	packedSize += utils::packedSize(stiffness);

	packedSize += 1;
	if (epointers != NULL) {
		packedSize += sizeof(size_t) + epointers->datatarray().size() * sizeof(int);
	}

	packedSize += utils::packedSize(eintervals);
	packedSize += utils::packedSize(eintervalsDistribution);

	packedSize += utils::packedSize(data.size());
	for (size_t i = 0; i < data.size(); i++) {
		packedSize += data[i]->packedSize();
	}

	return packedSize;
}

void ElementStore::packFull(char* &p) const
{
	utils::pack(distribution.threads, p);
	utils::pack(distribution.process, p);
	utils::pack(distribution.code, p);

	utils::pack(offset, p);
	utils::pack(inputOffset, p);
	utils::pack(nodes, p);
	utils::pack(centers, p);
	utils::pack(body, p);
	utils::pack(contact, p);
	utils::pack(material, p);
	utils::pack(regions, p);
	utils::pack(faceNeighbors, p);
	utils::pack(edgeNeighbors, p);
	utils::pack(volumeIndices, p);
	utils::pack(stiffness, p);

	utils::pack(epointers != NULL, p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());
		for (size_t i = 0; i < epointers->datatarray().size(); ++i) {
			eindices.push_back(static_cast<int>(epointers->datatarray()[i]->code));
		}
		utils::pack(eindices, p);
	}

	utils::pack(eintervals, p);
	utils::pack(eintervalsDistribution, p);

	utils::pack(data.size(), p);
	for (size_t i = 0; i < data.size(); i++) {
		data[i]->pack(p);
	}
}

void ElementStore::unpackFull(const char* &p)
{
	utils::unpack(distribution.threads, p);
	utils::unpack(distribution.process, p);
	utils::unpack(distribution.code, p);

	utils::unpack(offset, p);
	utils::unpack(inputOffset, p);
	utils::unpack(nodes, p);
	utils::unpack(centers, p);
	utils::unpack(body, p);
	utils::unpack(contact, p);
	utils::unpack(material, p);
	utils::unpack(regions, p);
	utils::unpack(faceNeighbors, p);
	utils::unpack(edgeNeighbors, p);
	utils::unpack(volumeIndices, p);
	utils::unpack(stiffness, p);

	bool notnull;
	utils::unpack(notnull, p);
	if (notnull) {
		std::vector<int> eindices;
		utils::unpack(eindices, p);
		if (epointers != NULL) {
			delete epointers;
		}
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(info::env::OMP_NUM_THREADS, distribution.process.size));
		for (esint i = 0; i < distribution.process.size; ++i) {
			epointers->datatarray()[i] = &Mesh::edata[eindices[i]];
		}
	}

	utils::unpack(eintervals, p);
	utils::unpack(eintervalsDistribution, p);

	size_t size;
	utils::unpack(size, p);
	for (size_t i = 0; i < size; i++) {
		data.push_back(new ElementData(p));
	}
}

size_t ElementStore::packedSize() const
{
	return
			nodes->packedSize() +
			sizeof(size_t) + epointers->datatarray().size() * sizeof(int) +
			utils::packedSize(body) +
			utils::packedSize(contact) +
			utils::packedSize(eintervals);
}

void ElementStore::pack(char* &p) const
{
	nodes->pack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());

		size_t threads = info::env::OMP_NUM_THREADS;
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution.threads[t]; i < distribution.threads[t + 1]; ++i) {
				eindices.push_back(static_cast<int>(epointers->datatarray()[i]->code));
			}
		}
		utils::pack(eindices, p);
	}
	utils::pack(body, p);
	utils::pack(contact, p);
	utils::pack(eintervals, p);
}

void ElementStore::unpack(const char* &p)
{
	if (nodes == NULL) {
		nodes = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
	}
	if (epointers == NULL) {
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, 0));
	}

	nodes->unpack(p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		utils::unpack(eindices, p);
		if (epointers != NULL) {
			delete epointers;
		}
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, distribution.process.size));
		for (esint i = 0; i < distribution.process.size; ++i) {
			epointers->datatarray()[i] = &Mesh::edata[eindices[i]];
		}
	}
	utils::unpack(body, p);
	utils::unpack(contact, p);
	utils::unpack(eintervals, p);
}

size_t ElementStore::packedDataHeaderSize() const
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

void ElementStore::packDataHeader(char* &p) const
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

void ElementStore::unpackDataHeader(const char* &p)
{
	size_t datasize;
	utils::unpack(datasize, p);
	for (size_t i = 0; i < datasize; i++) {
		data.push_back(new ElementData(0, NamedData::DataType::VECTOR, {}));
		utils::unpack(data[i]->dimension, p);
		utils::unpack(data[i]->dataType, p);
		utils::unpack(data[i]->name, p);
		utils::unpack(data[i]->restriction, p);
	}
}

size_t ElementStore::packedDataSize() const
{
	size_t size = 1;
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->name.size()) {
			size += utils::packedSize(data[i]->data);
		}
	}
	return size;
}

void ElementStore::packData(char* &p) const
{
	for (size_t i = 0; i < data.size(); i++) {
		if (data[i]->name.size()) {
			utils::pack(data[i]->data, p);
		}
	}
}

void ElementStore::unpackData(const char* &p)
{
	for (size_t i = 0; i < data.size(); i++) {
		utils::unpack(data[i]->data, p);
	}
}

ElementStore::~ElementStore()
{
	if (offset != NULL) { delete offset; }
	if (inputOffset != NULL) { delete inputOffset; }
	if (nodes != NULL) { delete nodes; }
	if (centers != NULL) { delete centers; }

	if (body != NULL) { delete body; }
	if (contact != NULL) { delete contact; }
	if (material != NULL) { delete material; }
	if (regions != NULL) { delete regions; }
	if (epointers != NULL) { delete epointers; }

	if (faceNeighbors != NULL) { delete faceNeighbors; }
	if (edgeNeighbors != NULL) { delete edgeNeighbors; }

	if (volumeIndices != NULL) { delete volumeIndices; }

	if (stiffness != NULL) { delete stiffness; }

	for (size_t i = 0; i < data.size(); i++) {
		delete data[i];
	}
}

void ElementStore::store(const std::string &file)
{
	std::ofstream os(file + std::to_string(info::mpi::rank) + ".txt");

	Store::storedata(os, "offset", offset);
	Store::storedata(os, "outputOffset", inputOffset);
	Store::storedata(os, "nodes", nodes);
	Store::storedata(os, "centers", centers);

	Store::storedata(os, "body", body);
	Store::storedata(os, "material", material);
	Store::storedata(os, "regions", regions);
	Store::storedata(os, "epointers", epointers);

	Store::storedata(os, "neighbors", faceNeighbors);
}

void ElementStore::permute(const std::vector<esint> &permutation, const std::vector<size_t> &threading)
{
	distribution.threads = threading;

	if (offset != NULL) { offset->permute(permutation, threading); }
	if (inputOffset != NULL) { inputOffset->permute(permutation, threading); }
	if (nodes != NULL) { nodes->permute(permutation, threading); }
	if (centers != NULL) { centers->permute(permutation, threading); }

	if (body != NULL) { body->permute(permutation, threading); }
	if (contact != NULL) { contact->permute(permutation, threading); }
	if (material != NULL) { material->permute(permutation, threading); }
	if (regions != NULL) { regions->permute(permutation, threading); }

	if (epointers != NULL) { epointers->permute(permutation, threading); }

	if (faceNeighbors != NULL) { faceNeighbors->permute(permutation, threading); }
	if (edgeNeighbors != NULL) { edgeNeighbors->permute(permutation, threading); }

	if (volumeIndices != NULL) { volumeIndices->permute(permutation, threading); }

	if (stiffness != NULL) { stiffness->permute(permutation, threading); }

	// TODO: permute data
}

void ElementStore::reindex(const serializededata<esint, esint> *nIDs)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = nodes->begin(t)->begin(); n != nodes->end(t)->begin(); ++n) {
			*n = std::lower_bound(nIDs->datatarray().begin(), nIDs->datatarray().end(), *n) - nIDs->datatarray().begin();
		}
	}
}

ElementData* ElementStore::appendData(int dimension, NamedData::DataType datatype, const std::string &name, step::TYPE restriction)
{
	this->data.push_back(new ElementData(dimension, datatype, name));
	data.back()->restriction = restriction;
	data.back()->data.resize(distribution.process.size * dimension);
	return this->data.back();
}

void ElementData::statistics(const tarray<esint> &elements, esint totalsize, Statistics *statistics) const
{
	for (int d = 0; d <= nstatistics(); d++) {
		(statistics + d)->reset();
	}

	esint doffset = nstatistics() != dimension ? 1 : 0;
	for (auto e = elements.begin(); e != elements.end(); ++e) {
		double value = 0;
		for (int d = 0; d < dimension; d++) {
			value += store[*e * dimension + d] * store[*e * dimension + d];
			(statistics + d + doffset)->min    = std::min((statistics + d + doffset)->min, store[*e * dimension + d]);
			(statistics + d + doffset)->max    = std::max((statistics + d + doffset)->max, store[*e * dimension + d]);
			(statistics + d + doffset)->avg   += store[*e * dimension + d];
			(statistics + d + doffset)->norm  += store[*e * dimension + d] * store[*e * dimension + d];
			(statistics + d + doffset)->absmin = std::min((statistics + d + doffset)->min, std::fabs(store[*e * dimension + d]));
			(statistics + d + doffset)->absmax = std::max((statistics + d + doffset)->max, std::fabs(store[*e * dimension + d]));
		}
		if (dataType == DataType::VECTOR) {
			value = std::sqrt(value);
			statistics->min = std::min(statistics->min, value);
			statistics->max = std::max(statistics->max, value);
			statistics->avg += value;
			statistics->norm += value * value;
			statistics->absmin = std::min(statistics->min, std::fabs(value));
			statistics->absmax = std::max(statistics->max, std::fabs(value));
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














