
#include "w.openvdb.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
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

void OpenVDBWrapper::add_grid(size_t distMax, size_t dataMax, esint *dist, _Point<short>* voxels, float *data, const std::string &name)
{
#ifdef HAVE_OPENVDB
	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
	grid->setGridClass(openvdb::GRID_LEVEL_SET);
	grid->setTransform(openvdb::math::Transform::createLinearTransform(_data->mat_));
	grid->setName(name);
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

	for (int r = 0; r < info::mpi::size; ++r) {
		for (esint e = 0; e < dist[0]; ++e) {
			for (esint v = dist[1 + e]; v < dist[1 + e + 1]; ++v) {
				openvdb::Coord xyz(voxels[v].x, voxels[v].y, voxels[v].z);
				accessor.setValue(xyz, data[e]);
			}
		}
		dist += distMax;
		voxels += dataMax;
		data += distMax;
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


