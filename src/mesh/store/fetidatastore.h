
#ifndef SRC_MESH_STORE_FETIDATASTORE_H_
#define SRC_MESH_STORE_FETIDATASTORE_H_

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct FETIDataStore {

	// B0 from kernels
	serializededata<esint, esint>* domainDual;

	// B0 from corners
	std::vector<esint> corners;

	// Regularization from fix points
	std::vector<esint> surfaceFixPoints, sFixPointsDistribution;
	std::vector<esint> innerFixPoints, iFixPointsDistribution;

	FETIDataStore();
	~FETIDataStore();

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);
};

}



#endif /* SRC_MESH_STORE_FETIDATASTORE_H_ */
