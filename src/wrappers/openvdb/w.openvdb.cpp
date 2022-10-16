
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
	openvdb::math::Mat4d mat;

	OpenVDBWrapperData(const Point &origin, const Point &size, const _Point<short> &density)
	{
		double x = size.x / density.x, y = size.y / density.y, z = size.z / density.z;
		mat = openvdb::math::Mat4d(
				x, 0.0, 0.0, 0.0,
				0.0, y, 0.0, 0.0,
				0.0, 0.0, z, 0.0,
				origin.x - .5 * x, origin.y - .5 * y, origin.z - .5 * z, 1.0);
	}

	~OpenVDBWrapperData()
	{

	}
};

struct OpenVDBFloatWrapper {
	OpenVDBFloatWrapper(openvdb::FloatGrid::Ptr grid): grid(grid) {}
	openvdb::FloatGrid::Ptr grid;
};

}
#endif


using namespace espreso;

OpenVDBWrapper::OpenVDBWrapper(const Point &origin, const Point &size, const _Point<short> &density)
{
#ifndef HAVE_OPENVDB
	eslog::warning("ESPRESO run-time warning: cannot store output to OpenVDB (the library is not linked).\n");
	_wrapper = nullptr;
#else
	_wrapper = new OpenVDBWrapperData(origin, size, density);
#endif
}

OpenVDBWrapper::~OpenVDBWrapper()
{
	for (size_t i = 0; i < _data.size(); ++i) {
		delete _data[i];
	}
#ifdef HAVE_OPENVDB
	delete _wrapper;
#endif
}

OpenVDBWrapper::FloatData* OpenVDBWrapper::addFloat(const std::string &name)
{
	FloatData *data = new FloatData(name);
	_data.push_back(data);
#ifdef HAVE_OPENVDB
	_wrapper->grids.push_back(data->wrapper->grid);
	_wrapper->grids.back()->setGridClass(openvdb::GRID_LEVEL_SET);
	_wrapper->grids.back()->setTransform(openvdb::math::Transform::createLinearTransform(_wrapper->mat));
	_wrapper->grids.back()->setName(name);
#endif
	return data;
}

OpenVDBWrapper::FloatData::FloatData(const std::string &name)
{
#ifdef HAVE_OPENVDB
	wrapper = new OpenVDBFloatWrapper(openvdb::FloatGrid::create());
#endif
}

OpenVDBWrapper::FloatData::~FloatData()
{
#ifdef HAVE_OPENVDB
	delete wrapper;
#endif
}

void OpenVDBWrapper::FloatData::insert(const VolumePacker &packer, char *voxels, float *values)
{
#ifdef HAVE_OPENVDB
	openvdb::FloatGrid::Accessor accessor = wrapper->grid->getAccessor();

	packer.unpack(voxels, [&] (const _Point<short> &voxel, const size_t &vindex) {
		openvdb::Coord xyz(voxel.x, voxel.y, voxel.z);
		accessor.setValue(xyz, values[vindex]);
	});
#endif
}

void OpenVDBWrapper::store(const std::string &file)
{
#ifdef HAVE_OPENVDB
	openvdb::io::File(file).write(_wrapper->grids);
	for (size_t i = 0; i < _wrapper->grids.size(); ++i) {
		_wrapper->grids[i]->clear();
	}
#endif
}


