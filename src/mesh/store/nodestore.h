
#ifndef SRC_MESH_STORE_NODESTORE_H_
#define SRC_MESH_STORE_NODESTORE_H_

#include "nameddata.h"
#include "nodeuniquenessinfo.h"
#include "basis/containers/point.h"

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

template <typename TEData> class tarray;
template <typename TEBoundaries, typename TEData> class serializededata;
struct Statistics;

struct NodeData: public NamedData {
	NodeData(int dimension, DataType datatype, const std::string &name): NamedData(dimension, datatype, name) {}
	NodeData(const char* &packedData): NamedData(packedData) {}

	void statistics(const tarray<esint> &nodes, esint totalsize, Statistics *statistics) const;
	void synchronize();
	static void synchronize(std::vector<double> &data, int dimension);
};

// store nodes held by lower ranks first, them all my nodes
// i.e. nodes can be divided into two intervals: <0, nhalo> <nhalo, size>; second sorted according IDs
struct NodeStore {
	friend class Mesh;

	void store(const std::string &file);

	void permute(const std::vector<esint> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution);

	NodeData* appendData(int dimension, NamedData::DataType datatype, const std::string &name = "", step::TYPE restriction = step::TYPE::TIME, bool toOutput = true);

	std::vector<esint> gatherNodeDistribution();
	std::vector<esint> gatherUniqueNodeDistribution();

	esint size;
	std::vector<size_t> distribution;

	serializededata<esint, esint>* IDs;
	serializededata<esint, esint>* elements;

	serializededata<esint, Point>* originCoordinates;
	serializededata<esint, Point>* coordinates;
	serializededata<esint, int>* ranks;
	serializededata<esint, int>* domains;

	NodeUniquenessInfo uniqInfo;

	std::vector<NodeData*> data;

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

	NodeStore();
	~NodeStore();
};

}


#endif /* SRC_MESH_STORE_NODESTORE_H_ */
