
#ifndef SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_
#define SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_

#include "basis/containers/point.h"

#include <string>
#include <vector>

namespace espreso {

struct OpenVDBWrapperData;
struct OpenVDBFloatWrapper;

struct OpenVDBWrapper {

	struct Data { virtual ~Data() {} };

	struct FloatData: Data {
		FloatData(const std::string &name);
		~FloatData();
		void insert(const _Point<short> &voxel, const float &value);

		OpenVDBFloatWrapper *wrapper = nullptr;
	};

	OpenVDBWrapper(const Point &origin, const Point &size, const _Point<short> &density);
	~OpenVDBWrapper();

	FloatData* addFloat(const std::string &name);
	void store(const std::string &file);

protected:
	OpenVDBWrapperData *_wrapper = nullptr;
	std::vector<Data*> _data;
};

}

#endif /* SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_ */
