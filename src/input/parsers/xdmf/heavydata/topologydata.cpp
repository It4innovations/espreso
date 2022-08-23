
#include "topologydata.h"
#include "basis/containers/tarray.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/element.h"
#include "input/meshbuilder.h"
#include "input/parsers/xdmf/lightdata/xdmftopology.h"
#include "input/parsers/xdmf/lightdata/xdmfdataitem.h"
#include "wrappers/hdf5/w.hdf5.h"

#include <numeric>

using namespace espreso;

int TopologyData::firstrank = 0;

static Element::CODE recognize(XDMFTopology::Type type)
{
	// some elements types are not supported by ESPRESO
	switch (type) {
	case XDMFTopology::Type::Polyvertex:      return Element::CODE::POINT1;
	case XDMFTopology::Type::Polyline:        return Element::CODE::LINE2;
	case XDMFTopology::Type::Triangle:        return Element::CODE::TRIANGLE3;
	case XDMFTopology::Type::Quadrilateral:   return Element::CODE::SQUARE4;
	case XDMFTopology::Type::Tetrahedron:     return Element::CODE::TETRA4;
	case XDMFTopology::Type::Pyramid:         return Element::CODE::PYRAMID5;
	case XDMFTopology::Type::Wedge:           return Element::CODE::PRISMA6;
	case XDMFTopology::Type::Hexahedron:      return Element::CODE::HEXA8;
	case XDMFTopology::Type::Edge_3:          return Element::CODE::LINE3;
	case XDMFTopology::Type::Triangle_6:      return Element::CODE::TRIANGLE6;
	case XDMFTopology::Type::Quadrilateral_8: return Element::CODE::SQUARE8;
	case XDMFTopology::Type::Tetrahedron_10:  return Element::CODE::TETRA10;
	case XDMFTopology::Type::Pyramid_13:      return Element::CODE::PYRAMID13;
	case XDMFTopology::Type::Wedge_15:        return Element::CODE::PRISMA15;
	case XDMFTopology::Type::Hexahedron_20:   return Element::CODE::HEXA20;
	default:
		return Element::CODE::SIZE;
	}
}

TopologyData::TopologyData(XDMFTopology *topology, XDMFDataItem *topologydata)
: etype((int)Element::CODE::SIZE), esize(1)
{
	switch (topology->type) {
	case XDMFTopology::Type::CoRectMesh2D:
	case XDMFTopology::Type::CoRectMesh3D:
	case XDMFTopology::Type::RectMesh2D:
	case XDMFTopology::Type::RectMesh3D:
	case XDMFTopology::Type::SMesh2D:
	case XDMFTopology::Type::SMesh3D:
		eslog::error("XDMF parser error: unsupported topology type.\n");
		break;
	default: break;
	}

	name = Parser::split(topologydata->data, ":")[1];
	etype = (int)recognize(topology->type);
	esize = topologydata->dimensions.size() > 1 ? topologydata->dimensions[1] : 1;
	if (topologydata->dimensions[0] * esize * sizeof(esint) < 1024 * 1024) {
		// in order to increase change to correct recognition of elements
		distribution.reserve(info::mpi::size + 1);
		for (int i = 0; i <= info::mpi::size; ++i) {
			if (i <= firstrank) {
				distribution.push_back(0);
			} else {
				distribution.push_back(topologydata->dimensions[0]);
			}
		}
		firstrank = (firstrank + 1) % info::mpi::size;
	} else {
		distribution = tarray<esint>::distribute(info::mpi::size, topologydata->dimensions[0], 5 * align);
	}
}

void TopologyData::read(HDF5 &file)
{
	if (MPITools::subset->across.rank == 0) {
		esint nelements = distribution[info::mpi::rank + MPITools::subset->across.size] - distribution[info::mpi::rank];
		data.reserve(esize * nelements + align);
		data.resize(esize * nelements);
		file.read(name.c_str(), HDF5::ESINT, data.data(), esize, nelements, distribution[info::mpi::rank]);
	}
}













