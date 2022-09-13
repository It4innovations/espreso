
#ifndef SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_
#define SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_

#include <vector>
#include <string>

namespace espreso {

struct OpenVDBWrapperData;

template <typename TType> class _Point;
template <typename TEBoundaries, typename TEData> class serializededata;
struct ElementData;

struct OpenVDBWrapper {

	OpenVDBWrapper();
	~OpenVDBWrapper();

	void add_grid(size_t distMax, size_t dataMax, esint *dist, _Point<short>* voxels, float *data, const std::string &name);
	void store_grids(const char *name);

protected:
	OpenVDBWrapperData *_data;
};

}

#endif /* SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_ */
