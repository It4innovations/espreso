
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

	OpenVDBWrapper(const std::string &path, const std::string &directory, const std::string &name);
	~OpenVDBWrapper();

	void add_grid(std::vector<int> voxel_indices_dist, std::vector<_Point<int>> voxel_indices, std::vector<double> data, std::string data_name, int data_dimension);
	void store_grids(const char *name);

protected:
	OpenVDBWrapperData *_data;
	std::string _path, _directory, _name;
};

}

#endif /* SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_ */
