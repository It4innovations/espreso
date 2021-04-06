
#ifndef SRC_BASIS_UTILITIES_PACKING_H_
#define SRC_BASIS_UTILITIES_PACKING_H_

#include "basis/containers/serializededata.h"
#include "math/vector.dense.h"

#include <string>
#include <cstring>
#include <vector>

namespace espreso {
namespace utils {

	template <typename TEBoundaries, typename TEData>
	size_t packedSize(serializededata<TEBoundaries, TEData> *data);

	size_t packedSize(const VectorDense &data);

	template<typename Ttype>
	size_t packedSize(const Ttype &data);

	template<typename Ttype>
	size_t packedSize(const std::vector<Ttype> &data);

	template <typename TEBoundaries, typename TEData>
	void pack(serializededata<TEBoundaries, TEData> *data, char* &p);

	void pack(VectorDense &data, char* &p);

	template<typename Ttype>
	void pack(const Ttype &data, char* &p);

	template<typename Ttype>
	void pack(const std::vector<Ttype> &data, char* &p);

	template <typename TEBoundaries, typename TEData>
	void unpack(serializededata<TEBoundaries, TEData> *&data, const char* &p);

	void unpack(VectorDense &data, const char* &p);

	template<typename Ttype>
	void unpack(Ttype &data, const char* &p);

	template<typename Ttype>
	void unpack(std::vector<Ttype> &data, const char* &p);
}
}


#include "packing.hpp"


#endif /* SRC_BASIS_UTILITIES_PACKING_H_ */
