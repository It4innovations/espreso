
#include "xdmf.h"
#include "writer/xdmfwritter.h"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/xml.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "wrappers/hdf5/w.hdf5.h"

namespace espreso {

/******************************************************************************
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="3">
  <Domain>                                           // XDMFData::domain
    <Grid GridType="Collection">                     // XDMFData::collection
      <Grid GridType="Tree">                         // XDMFData::tree
        <Grid Name="region->name">                   // XDMFData::region[index]
          <Geometry> reference to domain </Geometry> // XDMFData::geometry[index]
          <Topology> reference to domain </Topology> // XDMFData::topology[index]
        </Grid>
      </Grid> // Tree
      .
      .
      .
      <DataItem> h5::geometry </DataItem> // XDMF::Geometry
      <DataItem> h5::topology </DataItem> // XDMF::Topology
      .
      .
    </Grid> // Collection
  </Domain>
</Xdmf>
******************************************************************************/

struct XDMFData {
	esint iteration;
	XML::Element *domain;
	XML::Element *collection;
	XML::Element *tree;
	std::vector<XML::Element*> region, geometry, topoloty;
};

}

using namespace espreso;

static Element::CODE getcode(const std::vector<esint> &ecounters) {
	Element::CODE code = Element::CODE::POINT1;
	for (size_t i = 0; i < ecounters.size(); ++i) {
		if (ecounters[i]) {
			if (code != Element::CODE::POINT1) {
				return Element::CODE::SIZE;
			}
			code = (Element::CODE)(i);
		}
	}
	return code;
};

static esint getoffset(const std::vector<esint> &eoffsets, Element::CODE code) {
	esint offset = 0;
	for (size_t i = 0; i < eoffsets.size(); ++i) {
		offset += eoffsets[i] * Mesh::edata[i].nodes;
		if (code == Element::CODE::SIZE) {
			offset += eoffsets[i];
		}
		if (i == (size_t)Element::CODE::LINE2) {
			offset += eoffsets[i];
		}
	}
	return offset;
};

static esint gettotalsize(const std::vector<esint> &ecounters, Element::CODE code) {
	return getoffset(ecounters, code);
};

static void fillGeometry(const std::string &path, XDMFData *xml, XDMF::Geometry &heavydata, const Mesh &mesh, const RegionStore *region, size_t rindex)
{
	heavydata.coordinates.reserve(region->nodeInfo.size * mesh.dimension);
	for (esint n = 0, i = region->nodeInfo.nhalo; n < region->nodeInfo.size; ++n, ++i) {
		for (int d = 0; d < mesh.dimension; ++d) {
			heavydata.coordinates.push_back(mesh.nodes->coordinates->datatarray()[region->nodes->datatarray()[i]][d]);
		}
	}
	heavydata.name = region->name + "_COORDINATES";
	heavydata.dimension = mesh.dimension;
	heavydata.offset = region->nodeInfo.offset;
	heavydata.size = region->nodeInfo.size;
	heavydata.totalsize = region->nodeInfo.totalSize;

	auto geometry = xml->region[rindex]->element("Geometry");
	auto reference = geometry->element("DataItem");
	auto dataitem = xml->domain->element("DataItem");

	geometry->attribute("GeometryType", mesh.dimension == 2 ? "XY" : "XYZ");

	reference->attribute("Reference", "XML");
	reference->value = dataitem->ref(heavydata.name);

	dataitem->attribute("Name", heavydata.name);
	dataitem->attribute("DataType", "Float");
	dataitem->attribute("Format", "HDF");
	dataitem->attribute("Dimensions", std::to_string(region->nodeInfo.totalSize) + " " + std::to_string(mesh.dimension));
	dataitem->value = path + ".h5:" + heavydata.name;
}

static void fillTopology(const std::string &path, XDMFData *xml, XDMF::Topology &heavydata, const RegionStore *region, Element::CODE code, size_t rindex)
{
	// heavy data are already stored
	heavydata.name = region->name + "_TOPOLOGY";
	heavydata.dimension = code == Element::CODE::SIZE ? 1 : code == Element::CODE::LINE2 ? 3 : Mesh::edata[(int)code].nodes;
	heavydata.offset = getoffset(region->eoffsets, code) / heavydata.dimension;
	heavydata.size = heavydata.topology.size() / heavydata.dimension;
	heavydata.totalsize = code == Element::CODE::POINT1 ? region->nodeInfo.totalSize : gettotalsize(region->ecounters, code) / heavydata.dimension;

	auto topology = xml->region[rindex]->element("Topology");
	auto reference = topology->element("DataItem");
	auto dataitem = xml->domain->element("DataItem");

	topology->attribute("TopologyType", XDMFWritter::etype(code));
	topology->attribute("NumberOfElements", std::to_string(code == Element::CODE::POINT1 ? region->nodeInfo.totalSize : region->totalsize));

	reference->attribute("Reference", "XML");
	reference->value = dataitem->ref(heavydata.name);

	dataitem->attribute("Name", heavydata.name);
	dataitem->attribute("DataType", sizeof(esint) == 4 ? "Int": "Long");
	dataitem->attribute("Format", "HDF");
	if (heavydata.dimension != 1) {
		dataitem->attribute("Dimensions", std::to_string(heavydata.totalsize) + " " + std::to_string(heavydata.dimension));
	} else {
		dataitem->attribute("Dimensions", std::to_string(heavydata.totalsize));
	}
	dataitem->value = path + ".h5:" + heavydata.name;
}

static void fillAttribute(const std::string &path, XML::Element *xml, XDMF::Attribute &heavydata, const std::string &name, const std::string &type)
{
	auto attribute = xml->element("Attribute");
	auto dataitem = attribute->element("DataItem");

	attribute->attribute("Name", name);
	attribute->attribute("Center", type);

	if (heavydata.dimension > 1) {
		attribute->attribute("AttributeType", "Vector");
		dataitem->attribute("Dimensions", std::to_string(heavydata.totalsize) + " " + std::to_string(3));
	} else {
		attribute->attribute("AttributeType", "Scalar");
		dataitem->attribute("Dimensions", std::to_string(heavydata.totalsize));
	}
	dataitem->attribute("DataType", "Float");
	dataitem->attribute("Format", "HDF");
	dataitem->value = path + ".h5:" + heavydata.name;
}

static void fillGeometryAttribute(const std::string &path, XML::Element *xml, XDMF::Attribute &heavydata, const RegionStore *store, const NamedData *data, int iteration)
{
	heavydata.name = store->name + "_" + data->name + "_" + std::to_string(iteration);
	heavydata.dimension = data->dimension > 1 ? 3 : 1;
	heavydata.offset = store->nodeInfo.offset;
	heavydata.size = store->nodeInfo.size;
	heavydata.totalsize = store->nodeInfo.totalSize;
	heavydata.values.clear();
	heavydata.values.reserve((data->dimension > 1 ? 3 : 1) * store->nodeInfo.size);
	for (auto n = store->nodes->datatarray().cbegin() + store->nodeInfo.nhalo; n != store->nodes->datatarray().cend(); ++n) {
		for (int d = 0; d < data->dimension; ++d) {
			heavydata.values.push_back(data->data[*n * data->dimension + d]);
		}
		if (data->dimension == 2) {
			heavydata.values.push_back(0);
		}
	}

	fillAttribute(path, xml, heavydata, data->name, "Node");
}

static void fillTopologyAttribute(const std::string &path, XML::Element *xml, XDMF::Attribute &heavydata, const ElementsRegionStore *store, const NamedData *data, int iteration)
{
	heavydata.name = store->name + "_" + data->name + "_" + std::to_string(iteration);
	heavydata.dimension = data->dimension > 1 ? 3 : 1;
	heavydata.offset = store->offset;
	heavydata.size = store->size;
	heavydata.totalsize = store->totalsize;
	heavydata.values.clear();
	heavydata.values.reserve((data->dimension > 1 ? 3 : 1) * store->elements->structures());
	for (auto e = store->elements->datatarray().cbegin(); e != store->elements->datatarray().cend(); ++e) {
		for (int d = 0; d < data->dimension; ++d) {
			heavydata.values.push_back(data->data[*e * data->dimension + d]);

		}
		if (data->dimension == 2) {
			heavydata.values.push_back(0);
		}
	}

	fillAttribute(path, xml, heavydata, data->name, "Cell");
}

XDMF::XDMF(const Mesh &mesh)
: Visualization(mesh), _hdf5(NULL), _xml(NULL), _data(NULL)
{
	if (!HDF5::islinked()) {
		eslog::globalerror("ESPRESO run-time error: link parallel HDF5 library in order to store data in XDMF format.\n");
	}
	_xml = new XML();
	_data = new XDMFData();

	_xml->root.name = "Xdmf";
	_xml->root.attribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
	_xml->root.attribute("Version", "3");

	_data->iteration = 0;
	_data->domain = _xml->root.element("Domain");
	_data->collection = _data->domain->element("Grid");
	_data->collection->attribute("GridType", "Collection");
	_data->collection->attribute("CollectionType", "Temporal");
	_data->collection->attribute("Name", _name);
}

XDMF::~XDMF()
{
	if (_hdf5) { delete _hdf5; }
	if (_xml) { delete _xml; }
	if (_data) { delete _data; }
}

void XDMF::updateMesh()
{
	if (_measure) { eslog::startln("XDMF: STARTED", "XDMF"); }
	profiler::syncstart("store_xdmf");

	size_t rindex = 0, regions = _mesh.elementsRegions.size() + _mesh.boundaryRegions.size() - 2;
	std::vector<Geometry> geometries(regions);
	std::vector<Topology> topologies(regions);

	if (info::mpi::irank == 0) {
		_data->tree = _data->collection->element("Grid")->attribute("GridType","Tree");

		_data->region.resize(regions);
		_data->geometry.resize(regions);
		_data->topoloty.resize(regions);

		for (size_t r = 1; r < _mesh.elementsRegions.size(); ++r, ++rindex) {
			const ElementsRegionStore *region = _mesh.elementsRegions[r];
			_data->region[rindex] = _data->tree->element("Grid");
			_data->region[rindex]->attribute("Name", region->name);

			fillGeometry(_directory + _name, _data, geometries[rindex], _mesh, region, rindex);

			Element::CODE code = getcode(region->ecounters); // Element::CODE::SIZE is used for representation of Mixed topology
			esint prev = 0;
			auto element = _mesh.elements->procNodes->cbegin();
			for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); prev = *e++) {
				element += *e - prev;
				if (code == Element::CODE::SIZE) {
					topologies[rindex].topology.push_back(XDMFWritter::ecode(_mesh.elements->epointers->datatarray()[*e]->code));
				}
				if (_mesh.elements->epointers->datatarray()[*e]->code == Element::CODE::LINE2) {
					topologies[rindex].topology.push_back(1);
				}
				for (auto n = element->begin(); n != element->end(); ++n) {
					topologies[rindex].topology.push_back(region->getPosition(*n));
				}
			}

			fillTopology(_directory + _name, _data, topologies[rindex], region, code, rindex);
		}

		for (size_t r = 1; r < _mesh.boundaryRegions.size(); ++r, ++rindex) {
			const BoundaryRegionStore *region = _mesh.boundaryRegions[r];
			_data->region[rindex] = _data->tree->element("Grid");
			_data->region[rindex]->attribute("Name", region->name);

			fillGeometry(_directory + _name, _data, geometries[rindex], _mesh, region, rindex);

			Element::CODE code = getcode(region->ecounters);
			if (region->dimension) {
				auto epointer = region->epointers->datatarray().cbegin();
				for (auto e = region->procNodes->begin(); e != region->procNodes->end(); ++e, ++epointer) {
					if (code == Element::CODE::SIZE) {
						topologies[rindex].topology.push_back(XDMFWritter::ecode((*epointer)->code));
					}
					if ((*epointer)->code == Element::CODE::LINE2) {
						topologies[rindex].topology.push_back(1);
					}
					for (auto n = e->begin(); n != e->end(); ++n) {
						topologies[rindex].topology.push_back(region->getPosition(*n));
					}
				}
			} else {
				for (esint i = 0; i < region->nodeInfo.size; ++i) {
					topologies[rindex].topology.push_back(region->nodeInfo.offset + i);
				}
			}
			fillTopology(_directory + _name, _data, topologies[rindex], region, code, rindex);
		}
//		if (_withDecomposition) {
//			// TODO:
//		}
	}

	profiler::synccheckpoint("serialize");
	if (_measure) { eslog::checkpointln("XDMF: GEOMETRY SERIALIZED"); }

	if (info::mpi::grank == 0) {
		_xml->store(_path + _name + ".xmf");
	}
	profiler::synccheckpoint("lightdata");
	Communication::barrier(MPITools::asynchronous);
	if (_measure) { eslog::checkpointln("XDMF: LIGHTDATA STORED"); }

	if (_hdf5 == NULL && MPITools::subset->within.rank == 0) {
		_hdf5 = new HDF5((_path + _directory + _name).c_str(), MPITools::subset->across, HDF5::MODE::WRITE);
	}
	profiler::synccheckpoint("init_hdf5");
	Communication::barrier(MPITools::asynchronous);
	if (_measure) { eslog::checkpointln("HDF5: HDF5 INITIALIZED"); }

	if (MPITools::subset->withinsize > 1) { // gather data to writers
		if (sizeof(int) != sizeof(float)) {
			eslog::error("MESIO internal error: cannot gather data to writers.\n");
		}
		std::vector<int> sdata, rdata;
		if (MPITools::subset->within.rank) { // without root data
			size_t totalsize = 0;
			for (size_t i = 0; i < geometries.size(); ++i) {
				totalsize += 1 + geometries[i].dimension * geometries[i].size;
			}
			for (size_t i = 0; i < topologies.size(); ++i) {
				totalsize += 1 + topologies[i].topology.size();
			}
			sdata.reserve(totalsize);
			for (size_t i = 0; i < geometries.size(); ++i) {
				sdata.push_back(geometries[i].size);
				sdata.insert(sdata.end(), reinterpret_cast<int*>(geometries[i].coordinates.data()), reinterpret_cast<int*>(geometries[i].coordinates.data() + geometries[i].dimension * geometries[i].size));
			}
			for (size_t i = 0; i < topologies.size(); ++i) {
				sdata.push_back(topologies[i].size);
				sdata.insert(sdata.end(), topologies[i].topology.begin(), topologies[i].topology.end());
			}
		}
		Communication::gatherUnknownSize(sdata, rdata, &MPITools::subset->within);
		if (MPITools::subset->within.rank == 0) {
			size_t offset = 0;
			for (int r = 1; r < MPITools::subset->within.size; r++) {
				for (size_t i = 0; i < geometries.size(); ++i) {
					int size = rdata[offset++];
					geometries[i].size += size;
					geometries[i].coordinates.insert(geometries[i].coordinates.end(), reinterpret_cast<float*>(rdata.data() + offset), reinterpret_cast<float*>(rdata.data() + offset + geometries[i].dimension * size));
					offset += geometries[i].dimension * size;
				}
				for (size_t i = 0; i < topologies.size(); ++i) {
					int size = rdata[offset++];
					topologies[i].size += size;
					topologies[i].topology.insert(topologies[i].topology.end(), rdata.data() + offset, rdata.data() + offset + topologies[i].dimension * size);
					offset += topologies[i].dimension * size;
				}
			}
		}
	}

	profiler::synccheckpoint("reorder");
	if (_measure) { eslog::checkpointln("XDMF: GEOMETRY DATA REORDERED"); }

	if (MPITools::subset->within.rank == 0) {
		for (size_t i = 0; i < geometries.size(); ++i) {
			_hdf5->append(geometries[i].name.c_str(), HDF5::FLOAT, HDF5::FLOAT, geometries[i].coordinates.data(), geometries[i].dimension, geometries[i].size, geometries[i].offset, geometries[i].totalsize);
		}
		for (size_t i = 0; i < topologies.size(); ++i) {
			_hdf5->append(topologies[i].name.c_str(), HDF5::INT, HDF5::INT, topologies[i].topology.data(), topologies[i].dimension, topologies[i].size, topologies[i].offset, topologies[i].totalsize);
		}
	}
	profiler::synccheckpoint("write");
	Communication::barrier(MPITools::asynchronous);
	profiler::syncend("store_xdmf");
	if (_measure) { eslog::endln("XDMF: GEOMETRY STORED"); }
}

void XDMF::updateSolution()
{
	if (step::type == step::TYPE::FTT) {
		return; // TODO
	}

	size_t rindex = 0;
	if (_data->iteration) {
		_data->tree = _data->collection->element("Grid")->attribute("GridType","Tree");
		for (size_t r = 1; r < _mesh.elementsRegions.size(); ++r, ++rindex) {
			_data->region[rindex] = _data->tree->element("Grid");
			_data->region[rindex]->attribute("Name", _mesh.elementsRegions[r]->name);
			_data->geometry[rindex] = _data->region[rindex]->element(_data->geometry[rindex]);
			_data->topoloty[rindex] = _data->region[rindex]->element(_data->topoloty[rindex]);
		}
		for (size_t r = 1; r < _mesh.boundaryRegions.size(); ++r, ++rindex) {
			_data->region[rindex] = _data->tree->element("Grid");
			_data->region[rindex]->attribute("Name", _mesh.boundaryRegions[r]->name);
			_data->geometry[rindex] = _data->region[rindex]->element(_data->geometry[rindex]);
			_data->topoloty[rindex] = _data->region[rindex]->element(_data->topoloty[rindex]);
		}
	}
	_data->tree->element("Time")->attribute("Value", std::to_string(step::time::current));

	std::vector<XDMF::Attribute> attributes;
	attributes.reserve(_mesh.elements->data.size() + _mesh.nodes->data.size());
	for (size_t di = 0; di < _mesh.elements->data.size(); ++di) {
		if (storeData(_mesh.elements->data[di])) {
			rindex = 0;
			for (size_t r = 1; r < _mesh.elementsRegions.size(); ++r, ++rindex) {
				attributes.push_back({});
				fillTopologyAttribute(_directory + _name, _data->region[rindex], attributes.back(), _mesh.elementsRegions[r], _mesh.elements->data[di], _data->iteration);
			}
		}
	}

	for (size_t di = 0; di < _mesh.nodes->data.size(); ++di) {
		if (storeData(_mesh.nodes->data[di])) {
			rindex = 0;
			for (size_t r = 1; r < _mesh.elementsRegions.size(); ++r, ++rindex) {
				attributes.push_back({});
				fillGeometryAttribute(_directory + _name, _data->region[rindex], attributes.back(), _mesh.elementsRegions[r], _mesh.nodes->data[di], _data->iteration);
			}
			for (size_t r = 1; r < _mesh.boundaryRegions.size(); ++r, ++rindex) {
				attributes.push_back({});
				fillGeometryAttribute(_directory + _name, _data->region[rindex], attributes.back(), _mesh.boundaryRegions[r], _mesh.nodes->data[di], _data->iteration);
			}
		}
	}

	if (MPITools::subset->withinsize > 1) { // gather data to writers
		std::vector<float> sdata, rdata;
		if (MPITools::subset->within.rank) { // without root data
			size_t totalsize = 0;
			for (size_t i = 0; i < attributes.size(); ++i) {
				totalsize += 1 + attributes[i].values.size();
			}
			for (size_t i = 0; i < attributes.size(); ++i) {
				sdata.push_back(attributes[i].size);
				sdata.insert(sdata.end(), attributes[i].values.begin(), attributes[i].values.end());
			}
		}
		Communication::gatherUnknownSize(sdata, rdata, &MPITools::subset->within);
		if (MPITools::subset->within.rank == 0) {
			size_t offset = 0;
			for (size_t i = 0; i < attributes.size(); ++i) {
				esint size = rdata[offset++];
				attributes[i].size += size;
				attributes[i].values.insert(attributes[i].values.end(), rdata.data() + offset, rdata.data() + offset + attributes[i].dimension * size);
				offset += attributes[i].dimension * size;
			}
		}
	}

	if (MPITools::subset->within.rank == 0) {
		for (size_t i = 0; i < attributes.size(); ++i) {
			_hdf5->append(attributes[i].name.c_str(), HDF5::FLOAT, HDF5::FLOAT, attributes[i].values.data(), attributes[i].dimension, attributes[i].size, attributes[i].offset, attributes[i].totalsize);
		}
	}

	if (info::mpi::grank == 0) {
		_xml->store(_path + _name + ".xmf");
	}
	++_data->iteration;
}


