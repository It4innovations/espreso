
#ifndef SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_
#define SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_

#include "basis/containers/point.h"

#include <string>
#include <vector>

namespace espreso {

struct OpenVDBWrapperData;
struct OpenVDBFloatWrapper;

struct OpenVDBWrapper {

	struct Data {};

	struct FloatData: Data {
		FloatData(const std::string &name);
		void insert(esint elements, esint *dist, _Point<short>* voxels, float *data);

		OpenVDBFloatWrapper *wrapper;
	};

	OpenVDBWrapper(const Point &origin, const Point &size, const _Point<short> &density);
	~OpenVDBWrapper();

	FloatData* addFloat(const std::string &name);
	void store(const std::string &file);

protected:
	OpenVDBWrapperData *_wrapper;
	std::vector<Data*> _data;
};

}

#endif /* SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_ */
