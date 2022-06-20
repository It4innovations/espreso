
#include "w.openvdb.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "mesh/store/elementstore.h"

#ifdef HAVE_OPENVDB
#include <openvdb/openvdb.h>

namespace espreso {
struct OpenVDBWrapperData {
	openvdb::FloatGrid::Ptr grid;

	OpenVDBWrapperData()
	{
		grid = openvdb::FloatGrid::create();
		grid->setGridClass(openvdb::GRID_LEVEL_SET);
		grid->setName("density");
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

void OpenVDBWrapper::setTransformation()
{
#ifdef HAVE_OPENVDB
	// there should be a real transformation computed from the mesh
	openvdb::math::Mat4d mat_ = openvdb::math::Mat4d(
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0);

	_data->grid->setTransform(openvdb::math::Transform::createLinearTransform(mat_));
#endif
}

void OpenVDBWrapper::store(const char *name, const serializededata<esint, _Point<int> > *emap, const ElementData *data)
{
#ifdef HAVE_OPENVDB
	openvdb::FloatGrid::Accessor accessor = _data->grid->getAccessor();

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

	openvdb::io::File(name).write({ _data->grid });
#endif
}


