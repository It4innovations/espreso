
#include "w.openvdb.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "mesh/store/elementstore.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>

#ifdef HAVE_OPENVDB
#include <openvdb/openvdb.h>

namespace espreso {
struct OpenVDBWrapperData {
	openvdb::GridPtrVec grids;
	openvdb::math::Mat4d mat_;

	OpenVDBWrapperData()
	{
		mat_ = openvdb::math::Mat4d(
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0);
	}

	~OpenVDBWrapperData()
	{

	}
};
}
#endif


using namespace espreso;

OpenVDBWrapper::OpenVDBWrapper()
{
#ifndef HAVE_OPENVDB
	eslog::warning("ESPRESO run-time warning: cannot store output to OpenVDB (the library is not linked).\n");
	_data = nullptr;
#else
	_data = new OpenVDBWrapperData();
#endif
}

OpenVDBWrapper::~OpenVDBWrapper()
{
#ifdef HAVE_OPENVDB
	delete _data;
#endif
}

void OpenVDBWrapper::add_grid(std::vector<int> voxel_indices_dist, std::vector<_Point<int>> voxel_indices, std::vector<double> data, std::string data_name, int data_dimension)
{
#ifdef HAVE_OPENVDB
	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
	grid->setGridClass(openvdb::GRID_LEVEL_SET);
	grid->setTransform(openvdb::math::Transform::createLinearTransform(_data->mat_));
	grid->setName(data_name);
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

	auto value = data.cbegin();
	int max_voxel_inx = 0;
	int new_max_voxel_inx = 0;
	int false_voxels_count = 0;
	for (int voxels_dist_inx = 0; voxels_dist_inx < voxel_indices_dist.size() - 1; voxels_dist_inx++, value += data_dimension) {
		if(voxel_indices_dist[voxels_dist_inx] != -1){
			if(voxel_indices_dist[voxels_dist_inx] + max_voxel_inx > new_max_voxel_inx){
				new_max_voxel_inx = voxel_indices_dist[voxels_dist_inx] + max_voxel_inx;
			}
			if(voxel_indices_dist[voxels_dist_inx] == 0 && new_max_voxel_inx > max_voxel_inx){
				max_voxel_inx = new_max_voxel_inx;
			}

			int start_voxel_inx = voxel_indices_dist[voxels_dist_inx] + max_voxel_inx + false_voxels_count;
			int end_voxel_inx = voxel_indices_dist[voxels_dist_inx + 1] + max_voxel_inx + false_voxels_count;

			int true_inx = start_voxel_inx;
			for(int voxel_inx = start_voxel_inx; voxel_inx < end_voxel_inx; voxel_inx++, true_inx++){
				while(voxel_indices[true_inx].x == -1){
					false_voxels_count++;
					true_inx++;
				}
				openvdb::Coord xyz(voxel_indices[true_inx].x, voxel_indices[true_inx].y, voxel_indices[true_inx].z);
				accessor.setValue(xyz, (float)*value);
			}

			//	for (int d = 0; d < data_dimension; ++d) { // is it possible to store more dimensional data?
					// printf("%d %d %d -> %f\n", voxel->x, voxel->y, voxel->z, *value);
			//	}
			//}
		}		
	}

//	grid->pruneGrid(0.01); // increase sparseness
	_data->grids.push_back(grid);
#endif
}

void OpenVDBWrapper::store_grids(const char *name)
{
#ifdef HAVE_OPENVDB
	openvdb::io::File(name).write(_data->grids);
#endif
}


