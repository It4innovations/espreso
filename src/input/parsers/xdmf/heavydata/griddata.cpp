
#include "griddata.h"
#include "basis/utilities/parser.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"
#include "input/meshbuilder.h"
#include "input/parsers/mixedelementsparser.h"
#include "input/parsers/xdmf/lightdata/lightdata.h"
#include "input/parsers/xdmf/lightdata/xdmfdataitem.h"
#include "input/parsers/xdmf/lightdata/xdmfgrid.h"
#include "input/parsers/xdmf/lightdata/xdmfdomain.h"
#include "input/parsers/xdmf/lightdata/xdmfgeometry.h"
#include "input/parsers/xdmf/lightdata/xdmftopology.h"
#include "mesh/element.h"
#include "wrappers/hdf5/w.hdf5.h"

#include <vector>
#include <numeric>

using namespace espreso;

GridData::GridData(LightData &lightdata)
: _lightdata(lightdata)
{

}

void GridData::scan()
{
	_lightdata.recurse([&] (XDMFElement *e, XDMFElement::EType type) {
		switch (type) {
		case XDMFElement::EType::Grid:
			switch (dynamic_cast<XDMFGrid*>(e)->type) {
			case XDMFGrid::Type::Collection: break;
			case XDMFGrid::Type::Tree: break;
			case XDMFGrid::Type::Uniform: _grids.push_back(_Grid{ dynamic_cast<XDMFGrid*>(e), NULL, NULL, NULL, NULL }); break;
			case XDMFGrid::Type::Subset: break;
			}
			break;
		case XDMFElement::EType::Geometry: _grids.back().geometry = dynamic_cast<XDMFGeometry*>(e); break;
		case XDMFElement::EType::Topology: _grids.back().topology = dynamic_cast<XDMFTopology*>(e); break;
		default: break;
		}
	});

	auto xpath = [&] (XDMFDataItem *item) -> XDMFDataItem* {
		if (item->reference.size()) {
			if (StringCompare::caseInsensitiveEq(item->reference, "XML")) {
				std::vector<std::string> path = Parser::split(item->data, "[]");
				std::string name = Parser::split(path[1], "=")[1];
				name = name.substr(1, name.size() - 2);
				path = Parser::split(path[0], "/");

				if (
						path.size() == 4 &&
						StringCompare::caseInsensitiveEq(path[0], "") &&
						StringCompare::caseInsensitiveEq(path[1], "XDMF") &&
						StringCompare::caseInsensitiveEq(path[2], "Domain") &&
						StringCompare::caseInsensitiveEq(path[3], "DataItem")) {

					XDMFElement *e = _lightdata.domain.front()->get(name);
					if (dynamic_cast<XDMFDataItem*>(e)) {
						return dynamic_cast<XDMFDataItem*>(e);
					}
				}
				eslog::error("XDMF reader error: not supported DataItem reference.\n");
			} else {
				eslog::error("XDMF reader error: not supported reference type (only XML is supported).\n");
			}
		}
		return item;
	};

	_geometry.reserve(_grids.size());
	_topology.reserve(_grids.size());

	for (auto it = _grids.begin(); it != _grids.end(); ++it) {
		if (it->geometry->dataitem.size() != 1) {
			eslog::error("XDMF parser error: Geometry element with exactly one DataItem is supported.\n");
		}
		if (it->topology->dataitem.size() != 1) {
			eslog::error("XDMF parser error: Topology element with exactly one DataItem is supported.\n");
		}

		it->geometrydata = xpath(it->geometry->dataitem.front());
		it->topologydata = xpath(it->topology->dataitem.front());
		_geometry.push_back({ it->geometry, it->geometrydata });
		_topology.push_back({ it->topology, it->topologydata });
	}
}

void GridData::read()
{
	if (_grids.size() == 0) {
		return;
	}
	eslog::startln("HDF5: STARTED", "HDF5");

	std::string hdffile = _lightdata.dir + "/" + Parser::split(_grids.front().geometrydata->data, ":")[0];
	HDF5 hdf5(hdffile.c_str(), MPITools::subset->across, HDF5::MODE::READ);
	eslog::checkpointln("HDF5: INITIALIZED");

	for (size_t i = 0; i < _grids.size(); ++i) {
		_geometry[i].read(hdf5);
		_topology[i].read(hdf5);
	}
	eslog::checkpointln("HDF5: READ");

	int reduction = MPITools::subset->within.size;

	if (reduction > 1) { // scatter data to non-readers
		std::vector<size_t> displacement(reduction + 1);
		std::vector<char > sBuffer, rBuffer;

		int writer = info::mpi::rank - MPITools::subset->within.rank;

		size_t totalsize = 0;
		for (size_t i = 0; i < _grids.size(); ++i) {
			totalsize += _geometry[i].dimension * (_geometry[i].distribution[writer + reduction] - _geometry[i].distribution[writer]);
			totalsize += _topology[i].esize * ( _topology[i].distribution[writer + reduction] - _topology[i].distribution[writer]);
		}
		for (int r = writer; r < writer + reduction; ++r) {
			for (size_t i = 0; i < _grids.size(); ++i) {
				{
					size_t chunkoffset = _geometry[i].dimension * (_geometry[i].distribution[r] - _geometry[i].distribution[writer]);
					size_t chunksize = _geometry[i].dimension * (_geometry[i].distribution[r + 1] - _geometry[i].distribution[r]);
					displacement[r - writer + 1] += sizeof(float) * (chunkoffset + chunksize);
					if (MPITools::subset->within.rank == 0 && chunksize) {
						char* begin = reinterpret_cast<char*>(_geometry[i].data.data() + chunkoffset);
						char* end = reinterpret_cast<char*>(_geometry[i].data.data() + chunkoffset + chunksize);
						sBuffer.insert(sBuffer.end(), begin, end);
					}
				}
				{
					size_t chunkoffset = _topology[i].esize * ( _topology[i].distribution[r] - _topology[i].distribution[writer]);
					size_t chunksize = _topology[i].esize * (_topology[i].distribution[r + 1] - _topology[i].distribution[r]);
					displacement[r - writer + 1] += sizeof(esint) * (chunkoffset + chunksize);
					if (MPITools::subset->within.rank == 0 && chunksize) {
						char* begin = reinterpret_cast<char*>(_topology[i].data.data() + chunkoffset);
						char* end = reinterpret_cast<char*>(_topology[i].data.data() + chunkoffset + chunksize);
						sBuffer.insert(sBuffer.end(), begin, end);
					}
				}
			}
		}

		Communication::scatterv(sBuffer, rBuffer, displacement, &MPITools::subset->within);

		for (size_t i = 0, roffset = 0; i < _grids.size(); ++i) {
			{
				size_t chunksize = _geometry[i].dimension * (_geometry[i].distribution[info::mpi::rank + 1] - _geometry[i].distribution[info::mpi::rank]);
				float* data = reinterpret_cast<float*>(rBuffer.data() + roffset);
				_geometry[i].data.assign(data, data + chunksize);
				roffset += sizeof(float) * chunksize;
			}
			{
				size_t chunksize = _topology[i].esize * (_topology[i].distribution[info::mpi::rank + 1] - _topology[i].distribution[info::mpi::rank]);
				esint* data = reinterpret_cast<esint*>(rBuffer.data() + roffset);
				_topology[i].data.reserve(chunksize + TopologyData::align);
				_topology[i].data.assign(data, data + chunksize);
				roffset += sizeof(esint) * chunksize;
			}
		}
		eslog::checkpointln("HDF5: DATA SCATTERED");
	}

	Communication::barrier();
	eslog::checkpointln("HDF5: SYNCHRONIZED");

	{ // align data
		std::vector<esint, initless_allocator<esint> > sBuffer, rBuffer;
		sBuffer.reserve(_topology.size() * TopologyData::align);
		rBuffer.reserve(_topology.size() * TopologyData::align);
		for (size_t i = 0; i < _topology.size(); ++i) {
			sBuffer.insert(sBuffer.end(), _topology[i].data.begin(), _topology[i].data.begin() + TopologyData::align);
		}
		Communication::receiveUpper(sBuffer, rBuffer);
		for (size_t i = 0; i < _topology.size(); ++i) {
			if (info::mpi::rank + 1 != info::mpi::size) {
				_topology[i].data.insert(_topology[i].data.end(), rBuffer.begin() + i * TopologyData::align, rBuffer.begin() + (i + 1) * TopologyData::align);
			} else {
				_topology[i].data.resize(_topology[i].data.size() + TopologyData::align, 0);
			}
		}
	}

	eslog::endln("HDF5: FINISHED");
}

static Element::CODE recognize(int id)
{
	// some elements types are not supported by ESPRESO
	switch (id) {
		case  1: return Element::CODE::POINT1;
		case  2: return Element::CODE::LINE2;
//		case  3: return Element::CODE::HEXA8;
		case  4: return Element::CODE::TRIANGLE3;
		case  5: return Element::CODE::SQUARE4;
		case  6: return Element::CODE::TETRA4;
		case  7: return Element::CODE::PYRAMID5;
		case  8: return Element::CODE::PRISMA6;
		case  9: return Element::CODE::HEXA8;
//		case 16: return Element::CODE::HEXA8;
		case 34: return Element::CODE::LINE3;
//		case 35: return Element::CODE::HEXA8;
		case 36: return Element::CODE::TRIANGLE6;
		case 37: return Element::CODE::SQUARE8;
		case 38: return Element::CODE::TETRA10;
		case 39: return Element::CODE::PYRAMID13;
		case 40: return Element::CODE::PRISMA15;
//		case 41: return Element::CODE::HEXA8;
		case 48: return Element::CODE::HEXA20;
//		case 49: return Element::CODE::HEXA8;
//		case 50: return Element::CODE::HEXA8;
		default:
			return Element::CODE::SIZE;
		}
}

void GridData::parse(MeshBuilder &mesh)
{
	auto enodes = [] (size_t index, esint id) {
		switch (id) {
		case  1: return  2;
		case  2: return  3;
//		case  3: return  0;
		case  4: return  3;
		case  5: return  4;
		case  6: return  4;
		case  7: return  5;
		case  8: return  6;
		case  9: return  8;
//		case 16: return  8;
		case 34: return  3;
//		case 35: return  8;
		case 36: return  6;
		case 37: return  8;
		case 38: return 10;
		case 39: return 13;
		case 40: return 15;
//		case 41: return  8;
		case 48: return 20;
//		case 49: return  8;
//		case 50: return  8;
		}
		return 0;
	};

	std::vector<int> d(_topology.size()), dim(_topology.size());
	for (size_t i = 0; i < _topology.size(); ++i) {
		if (_topology[i].distribution[info::mpi::rank] == 0 && _topology[i].distribution[info::mpi::rank + 1] != 0) {
			d[i] = _topology[i].data.front();
		}
	}
	Communication::allReduce(d.data(), dim.data(), _topology.size(), MPI_INT, MPI_MAX);

	MixedElementsParser mixedparser;
	for (size_t i = 0; i < _topology.size(); ++i) {
		if (_topology[i].etype == (int)Element::CODE::SIZE) {
			mixedparser.add(_topology[i].data.data(), _topology[i].data.size() - TopologyData::align);
		}
	}
	mixedparser.parse(enodes);

	for (size_t i = 0, j = 0; i < _topology.size(); ++i) {
		if (_topology[i].etype == (int)Element::CODE::SIZE) {
			if (mixedparser.invalid[j++] != info::mpi::size) {
				eslog::warning("XDMF parser: synchronization of region '%s'.\n", _grids[i].grid->name.c_str());
			}
		}
	}

	size_t csize = 0;
	for (size_t i = 0; i < _geometry.size(); ++i) {
		csize += _geometry[i].distribution[info::mpi::rank + 1] - _geometry[i].distribution[info::mpi::rank];
	}

	mesh.coordinates.reserve(csize);
	for (size_t i = 0, nodes = 0; i < _geometry.size(); ++i) {
		esint ncoordinates = _geometry[i].distribution[info::mpi::rank + 1] - _geometry[i].distribution[info::mpi::rank];
		if (_geometry[i].dimension == 3) {
			for (esint n = 0; n < ncoordinates; ++n) {
				mesh.coordinates.push_back(Point(_geometry[i].data[3 * n + 0], _geometry[i].data[3 * n + 1], _geometry[i].data[3 * n + 2]));
			}
		} else {
			for (esint n = 0; n < ncoordinates; ++n) {
				mesh.coordinates.push_back(Point(_geometry[i].data[2 * n + 0], _geometry[i].data[2 * n + 1], 0));
			}
		}

		for (esint n = 0; n < ncoordinates; ++n) {
			mesh.nIDs.push_back(nodes + _geometry[i].distribution[info::mpi::rank] + n);
		}
		nodes += _geometry[i].distribution.back();
	}

	size_t esize = 0;
	for (size_t i = 0; i < _topology.size(); ++i) {
		if (_topology[i].etype == (int)Element::CODE::SIZE) {
			esize += _topology[i].data.size() / 2;
		} else {
			if (_topology[i].esize != 1) {
				esize += _topology[i].distribution[info::mpi::rank + 1] - _topology[i].distribution[info::mpi::rank];
			}
		}
	}

	mesh.esize.reserve(esize);
	mesh.etype.reserve(esize);
	mesh.eIDs.reserve(esize);
	mesh.enodes.reserve(esize);
	for (size_t i = 0, j = 0, nodes = 0, elements = 0; i < _topology.size(); ++i) {
		if (_topology[i].etype == (int)Element::CODE::SIZE) {
			std::vector<esint> ids;
			for (size_t n = mixedparser.first[j]; n + TopologyData::align < _topology[i].data.size(); ++n) {
				Element::CODE code = recognize(_topology[i].data[n]);
				if (code == Element::CODE::POINT1) {
					ids.push_back(nodes + _topology[i].data[n + 2]);
					n += 2;
				} else if (code == Element::CODE::LINE2) { // XDMF has POLYLINE (we assume pattern 2 2 id id, ...)
					mesh.etype.push_back((int)Element::CODE::LINE2);
					mesh.esize.push_back(2);
					++n;
					for (int nn = 0; nn < 2; ++nn, ++n) {
						mesh.enodes.push_back(_topology[i].data[n + 1] + nodes);
					}
					ids.push_back(elements + mixedparser.offset[j]++);
				} else {
					mesh.etype.push_back((int)code);
					mesh.esize.push_back(enodes(i, _topology[i].data[n]));
					for (int nn = 0; nn < mesh.esize.back(); ++nn, ++n) {
						mesh.enodes.push_back(_topology[i].data[n + 1] + nodes);
					}
					ids.push_back(elements + mixedparser.offset[j]++);
				}
			}

			if (dim[i] == 1) {
				mesh.nregions[_grids[i].grid->name] = ids;
			} else {
				mesh.eIDs.insert(mesh.eIDs.end(), ids.begin(), ids.end());
				mesh.eregions[_grids[i].grid->name].swap(ids);
				elements += mixedparser.nelements[j];
			}
			++j;
		} else {
			if (_topology[i].esize == 1) {
				std::vector<esint> ids;
				ids.reserve(_topology[i].data.size() - TopologyData::align);
				for (size_t n = 0; n < _topology[i].data.size() - TopologyData::align; ++n) {
					ids.push_back(nodes + _topology[i].data[n]);
				}
				mesh.nregions[_grids[i].grid->name].swap(ids);
			} else {
				size_t start = mesh.enodes.size();
				size_t offset = _topology[i].distribution[info::mpi::rank];
				size_t size = _topology[i].distribution[info::mpi::rank + 1] - _topology[i].distribution[info::mpi::rank];
				mesh.eregions[_grids[i].grid->name].resize(size);
				std::iota(mesh.eregions[_grids[i].grid->name].begin(), mesh.eregions[_grids[i].grid->name].end(), elements + offset);
				mesh.eIDs.insert(mesh.eIDs.end(), mesh.eregions[_grids[i].grid->name].begin(), mesh.eregions[_grids[i].grid->name].end());
				mesh.etype.resize(mesh.etype.size() + size, _topology[i].etype);

				if (_topology[i].etype == (int)Element::CODE::LINE2) {
					mesh.esize.resize(mesh.esize.size() + size, 2);
					for (size_t n = 0; n + TopologyData::align < _topology[i].data.size();) {
						n++; // nnodes
						mesh.enodes.push_back(_topology[i].data[n++]);
						mesh.enodes.push_back(_topology[i].data[n++]);
					}
				} else {
					mesh.esize.resize(mesh.esize.size() + size, _topology[i].esize);
					mesh.enodes.insert(mesh.enodes.end(), _topology[i].data.begin(), _topology[i].data.end() - TopologyData::align);
				}
				elements += _topology[i].distribution.back();
				for (size_t n = start; n < mesh.enodes.size(); ++n) {
					mesh.enodes[n] += nodes;
				}
			}
		}
		nodes += _geometry[i].distribution.back();
	}
}

