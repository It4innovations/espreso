
#include "geometrydata.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "input/parsers/xdmf/lightdata/xdmfgeometry.h"
#include "input/parsers/xdmf/lightdata/xdmfdataitem.h"
#include "wrappers/hdf5/w.hdf5.h"

#include <vector>
#include <numeric>

using namespace espreso;

GeometryData::GeometryData(XDMFGeometry *geometry, XDMFDataItem *geometrydata)
: dimension(0)
{
	switch (geometry->type) {
	case XDMFGeometry::Type::XY:
	case XDMFGeometry::Type::XYZ:
		break;
	default:
		eslog::error("XDMF parser: not supported Geometry format.'\n");
	}

	dimension = geometrydata->dimensions[1];
	name = Parser::split(geometrydata->data, ":")[1];
	distribution = tarray<esint>::distribute(info::mpi::size, geometrydata->dimensions[0]);
}

void GeometryData::read(HDF5 &file)
{
	if (MPITools::subset->within.rank == 0) {
		esint ncoordinates = distribution[info::mpi::rank + MPITools::subset->within.size] - distribution[info::mpi::rank];
		data.resize(ncoordinates * dimension);
		file.read(name.c_str(), HDF5::FLOAT, data.data(), dimension, ncoordinates, distribution[info::mpi::rank]);
	}
}


