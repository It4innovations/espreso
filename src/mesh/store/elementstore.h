
#ifndef SRC_MESH_STORE_ELEMENTSTORE_H_
#define SRC_MESH_STORE_ELEMENTSTORE_H_

#include "basis/containers/point.h"
#include "info.h"
#include "elementsinterval.h"
#include "contactinfo.h"
#include "nameddata.h"

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

template <typename TEData> class tarray;
template <typename TEBoundaries, typename TEData> class serializededata;
struct Statistics;

struct ElementData: public NamedData {
	ElementData(int dimension, DataType datatype, const std::string &name): NamedData(dimension, datatype, name) {}
	ElementData(const char* &packedData): NamedData(packedData) {}

	void statistics(const tarray<esint> &elements, esint totalsize, Statistics *statistics) const;
};

struct ElementStore: UniqueDataInfo {

	void store(const std::string &file);

	void permute(const std::vector<esint> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution);

	void reindex(const serializededata<esint, esint> *nIDs);

	ElementData* appendData(int dimension, NamedData::DataType datatype, const std::string &name = "", step::TYPE restriction = step::TYPE::TIME);

	std::vector<size_t> distribution;

	serializededata<esint, esint>* IDs;
	serializededata<esint, esint>* nodes;
	serializededata<esint, Point>* centers;

	serializededata<esint, int>* body;
	serializededata<esint, ContactInfo>* contact;
	serializededata<esint, int>* material;
	serializededata<esint, esint>* regions;
	serializededata<esint, Element*>* epointers;

	serializededata<esint, esint>* faceNeighbors;
	serializededata<esint, esint>* edgeNeighbors;

	serializededata<esint, double>* stiffness;

	UniqueDataInfo bodies;
	esint nclusters;

	int regionMaskSize;
	std::vector<esint> ecounters;
	std::vector<ElementsInterval> eintervals;
	std::vector<esint> eintervalsDistribution;

	std::vector<ElementData*> data;

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	size_t packedDataHeaderSize() const;
	void packDataHeader(char* &p) const;
	void unpackDataHeader(const char* &p);

	size_t packedDataSize() const;
	void packData(char* &p) const;
	void unpackData(const char* &p);

	ElementStore();
	~ElementStore();
};

}

#endif /* SRC_MESH_STORE_ELEMENTSTORE_H_ */
