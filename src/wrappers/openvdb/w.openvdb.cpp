
#include "w.openvdb.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "mesh/store/elementstore.h"

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
	eslog::globalerror("ESPRESO run-time error: cannot store output to OpenVDB (the library is not linked).\n");
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

void OpenVDBWrapper::add_grid(const serializededata<esint, _Point<int> > *emap, const ElementData *data)
{
#ifdef HAVE_OPENVDB
	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
	grid->setGridClass(openvdb::GRID_LEVEL_SET);
	grid->setTransform(openvdb::math::Transform::createLinearTransform(_data->mat_));
	grid->setName(data->name);
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

	auto value = data->data.cbegin();
	for (auto voxels = emap->begin(); voxels != emap->end(); ++voxels, value += data->dimension) {
		for (auto voxel = voxels->begin(); voxel != voxels->end(); ++voxel) {
			openvdb::Coord xyz(voxel->x, voxel->y, voxel->z);
			accessor.setValue(xyz, (float)*value);

			for (int d = 0; d < data->dimension; ++d) { // is it possible to store more dimensional data?
				// printf("%d %d %d -> %f\n", voxel->x, voxel->y, voxel->z, *value);
			}
		}
	}

	_data->grids.push_back(grid);
#endif
}

void OpenVDBWrapper::store_grids(const char *name)
{
#ifdef HAVE_OPENVDB
	openvdb::io::File(name).write(_data->grids);
#endif
}


