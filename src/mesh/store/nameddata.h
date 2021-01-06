
#ifndef SRC_MESH_STORE_NAMEDDATA_H_
#define SRC_MESH_STORE_NAMEDDATA_H_

#include "esinfo/stepinfo.h"

#include <string>
#include <vector>

namespace espreso {

struct NamedData {
	static std::vector<std::string> coordinateSuffixes;
	static std::vector<std::string> numberSuffixes;

	enum class DataType {
		SCALAR,
		NUMBERED,
		VECTOR,
		TENSOR_SYMM,
		TENSOR_ASYM
	};

	int dimension;
	DataType dataType;
	std::string name;
	step::TYPE restriction;

	std::vector<double> data, buffer, &store;

	NamedData(int dimension, DataType datatype, const std::string &name);
	NamedData(const char* &packedData);

	bool withSuffixes() const;
	bool onlySuffixed() const;
	int nstatistics() const;

	void toBuffer();

	std::string suffix(int index) const;

	size_t packedSize();
	void pack(char *&p);
};

}



#endif /* SRC_MESH_STORE_NAMEDDATA_H_ */
