
#include "vtklegacy.h"
#include "writer/vtkwritter.h"
#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"
#include "mesh/store/contactinterfacestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"

#include <functional>

using namespace espreso;

VTKLegacy::VTKLegacy()
{
	_suffix = ".vtk";
}

VTKLegacy::~VTKLegacy()
{

}

void VTKLegacy::updateMesh()
{
	if (_measure) { eslog::startln("VTK LEGACY: STARTED", "VTK LEGACY"); }
	profiler::syncstart("store_vtk");

	_points.resize((size_t)info::mesh->nodes->size * 13 * 3 + 1);
	_esize.resize(2 * 20 + 1); // max is the largest element size
	_ecode.resize(3 * (int)Element::CODE::SIZE + 1);
	_cells.resize(info::mesh->elementsRegions.size() + info::mesh->boundaryRegions.size() + info::mesh->contactInterfaces.size() - 2);
	for (esint n = 0; n < info::mesh->nodes->size; ++n) {
		const Point &p = info::mesh->nodes->coordinates->datatarray()[n];
		sprintf(_points.data() + (size_t)n * 13 * 3, "%12.5e %12.5e %12.5e\n", p.x, p.y, p.z);
	}
	for (int n = 0; n < 20; ++n) {
		sprintf(_esize.data() + n * 2, "%d", n + 1);
	}
	for (int n = 0; n < (int)Element::CODE::SIZE; ++n) {
		sprintf(_ecode.data() + n * 3, "%2d\n", _writer.ecode(Mesh::edata[n].code));
	}

	int index = 0, intsize = 11;
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r, ++index) {
		esint nnodes = info::mesh->elementsRegions[r]->nodes->datatarray().size();
		if (info::mesh->elements->distribution.process.size < nnodes * 5) {
			_cells[index].resize(info::mesh->nodes->size * (intsize + sizeof(int)) + 1);
			for (esint n = 0; n < nnodes; ++n) {
				esint nn = info::mesh->elementsRegions[r]->nodes->datatarray()[n];
				int chars = sprintf(_cells[index].data() + nn * (intsize + sizeof(int)), " %d", (int)info::mesh->elementsRegions[r]->nodeInfo.position[n]);
				memcpy(_cells[index].data() + nn * (intsize + sizeof(int)) + intsize, &chars, sizeof(int));
			}
		} else {
			_cells[index].resize(nnodes * (intsize + sizeof(int)) + 1);
			for (esint n = 0; n < nnodes; ++n) {
				int chars = sprintf(_cells[index].data() + n * (intsize + sizeof(int)), " %d", (int)info::mesh->elementsRegions[r]->nodeInfo.position[n]);
				memcpy(_cells[index].data() + n * (intsize + sizeof(int)) + intsize, &chars, sizeof(int));
			}
		}
	}
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r, ++index) {
		esint nnodes = info::mesh->boundaryRegions[r]->nodes->datatarray().size();
		_cells[index].resize(nnodes * (intsize + sizeof(int)) + 1);
		for (esint n = 0; n < nnodes; ++n) {
			int chars = sprintf(_cells[index].data() + n * (intsize + sizeof(int)), " %d", (int)info::mesh->boundaryRegions[r]->nodeInfo.position[n]);
			memcpy(_cells[index].data() + n * (intsize + sizeof(int)) + intsize, &chars, sizeof(int));
		}
	}
	for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r, ++index) {
		esint nnodes = info::mesh->contactInterfaces[r]->nodes->datatarray().size();
		_cells[index].resize(nnodes * (intsize + sizeof(int)) + 1);
		for (esint n = 0; n < nnodes; ++n) {
			int chars = sprintf(_cells[index].data() + n * (intsize + sizeof(int)), " %d", (int)info::mesh->contactInterfaces[r]->nodeInfo.position[n]);
			memcpy(_cells[index].data() + n * (intsize + sizeof(int)) + intsize, &chars, sizeof(int));
		}
	}
	profiler::synccheckpoint("preserialize");

	index = 0;
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r, ++index) {
		insertHeader();
		insertPoints(info::mesh->elementsRegions[r]);
		esint ncells = insertElements(info::mesh->elementsRegions[r], _cells[index]);
		if (info::ecf->output.store_decomposition) {
			if (isRoot()) {
				_writer.celldata(ncells);
			}
			insertDecomposition(info::mesh->elementsRegions[r]);
		}
		_writer.commitFile(_path + _name + "." + info::mesh->elementsRegions[r]->name + _suffix);
	}

	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r, ++index) {
		insertHeader();
		insertPoints(info::mesh->boundaryRegions[r]);
		insertElements(info::mesh->boundaryRegions[r], _cells[index]);
		_writer.commitFile(_path + _name + "." + info::mesh->boundaryRegions[r]->name + _suffix);
	}
	for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r, ++index) {
		insertHeader();
		insertPoints(info::mesh->contactInterfaces[r]);
		insertElements(info::mesh->contactInterfaces[r], _cells[index]);
		_writer.commitFile(_path + _name + "." + info::mesh->contactInterfaces[r]->name + _suffix);
	}
	profiler::synccheckpoint("serialize");
	if (_measure) { eslog::checkpointln("VTK LEGACY: GEOMETRY SERIALIZED"); }

	_writer.reorder();
	profiler::synccheckpoint("reorder");
	if (_measure) { eslog::checkpointln("VTK LEGACY: GEOMETRY DATA REORDERED"); }

	_writer.write();
	profiler::synccheckpoint("write");
	Communication::barrier(MPITools::asynchronous);
	profiler::syncend("store_vtk");
	if (_measure) { eslog::endln("VTK LEGACY: GEOMETRY STORED"); }
}

void VTKLegacy::updateMonitors(const step::Step &step)
{
}

void VTKLegacy::updateSolution(const step::Step &step, const step::Time &time)
{
	updateSolution(step, _path + _directory, _name, "." + std::to_string(time.current) + _suffix);
}

void VTKLegacy::updateSolution(const step::Step &step, const step::Frequency &frequency)
{
	updateSolution(step, _path + _directory, _name, "." + std::to_string(frequency.current) + _suffix);
}

void VTKLegacy::updateSolution(const step::Step &step, const std::string &dir, const std::string &name, const std::string &suffix)
{
	int index = 0;
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r, ++index) {
		insertHeader();
		insertPoints(info::mesh->elementsRegions[r]);
		esint ncells = insertElements(info::mesh->elementsRegions[r], _cells[index]);
		if (isRoot()) {
			_writer.celldata(ncells);
		}
		if (info::ecf->output.store_decomposition) {
			insertDecomposition(info::mesh->elementsRegions[r]);
		}
		for (size_t i = 0; i < info::mesh->elements->data.size(); ++i) {
			insertData(
					step,
					info::mesh->elements->data[i],
					info::mesh->elementsRegions[r]->elements->datatarray().size(),
					info::mesh->elementsRegions[r]->elements->datatarray().data());
		}

		if (isRoot()) {
			_writer.pointdata(info::mesh->elementsRegions[r]->nodeInfo.totalSize);
		}
		for (size_t i = 0; i < info::mesh->nodes->data.size(); ++i) {
			insertData(
					step,
					info::mesh->nodes->data[i],
					info::mesh->elementsRegions[r]->nodeInfo.size,
					info::mesh->elementsRegions[r]->nodes->datatarray().data() + info::mesh->elementsRegions[r]->nodeInfo.nhalo);
		}
		_writer.commitFile(dir + name + "." + info::mesh->elementsRegions[r]->name + suffix);
	}

	auto boundary = [&] (const BoundaryRegionStore *region) {
		insertHeader();
		insertPoints(region);
		insertElements(region, _cells[index]);
		if (isRoot()) {
			_writer.pointdata(region->nodeInfo.totalSize);
		}
		for (size_t i = 0; i < info::mesh->nodes->data.size(); ++i) {
			insertData(step, info::mesh->nodes->data[i], region->nodeInfo.size, region->nodes->datatarray().data() + region->nodeInfo.nhalo);
		}
		_writer.commitFile(dir + name + "." + region->name + suffix);
	};

	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r, ++index) {
		boundary(info::mesh->boundaryRegions[r]);
	}
	for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r, ++index) {
		boundary(info::mesh->contactInterfaces[r]);
	}

	_writer.reorder();
	_writer.write();
}

void VTKLegacy::insertHeader()
{
	if (isRoot()) {
		_writer.description("# vtk DataFile Version 2.0\n");
		_writer.description("EXAMPLE\n");
		_writer.description("ASCII\n");
		_writer.description("DATASET UNSTRUCTURED_GRID\n");
	}
}

void VTKLegacy::insertPoints(const RegionStore *store)
{
	if (isRoot()) {
		_writer.points(store->nodeInfo.totalSize);
	}
	for (esint n = 0, i = store->nodeInfo.nhalo; n < store->nodeInfo.size; ++n, ++i) {
		_writer.insert(13 * 3, _points.data() + (size_t)store->nodes->datatarray()[i] * 13 * 3);
	}
	_writer.groupData();
}

esint VTKLegacy::insertElements(const ElementsRegionStore *store, const std::vector<char, initless_allocator<char> > &data)
{
	esint enodes = 0;
	for (size_t i = 0; i < store->distribution.code.size(); i++) {
		if (store->distribution.code[i].totalSize) {
			enodes += store->distribution.code[i].totalSize * Mesh::edata[i].nodes;
		}
	}

	if (isRoot()) {
		_writer.cells(store->distribution.process.totalSize, store->distribution.process.totalSize + enodes);
	}
	int intsize = 11;
	for (auto e = store->elements->datatarray().cbegin(); e != store->elements->datatarray().cend(); ++e) {
		esint nnodes = store->nodes->datatarray().size();
		auto element = info::mesh->elements->nodes->cbegin() + *e;
		_writer.insert(element->size() > 9 ? 2 : 1, _esize.data() + element->size() * 2 - 2);
		if (info::mesh->elements->distribution.process.size < nnodes * 5) {
			for (auto n = element->begin(); n != element->end(); ++n) {
				int size = *reinterpret_cast<const int*>(data.data() + (intsize + sizeof(int)) * *n + intsize);
				_writer.insert(size, data.data() + (intsize + sizeof(int)) * *n);
			}
		} else {
			for (auto n = element->begin(); n != element->end(); ++n) {
				esint p = std::lower_bound(store->nodes->datatarray().begin(), store->nodes->datatarray().end(), *n) - store->nodes->datatarray().begin();
				int size = *reinterpret_cast<const int*>(data.data() + (intsize + sizeof(int)) * p + intsize);
				_writer.insert(size, data.data() + (intsize + sizeof(int)) * p);
			}
		}
		_writer.push('\n');
	}
	_writer.groupData();

	if (isRoot()) {
		_writer.celltypes(store->distribution.process.totalSize);
	}
	for (auto e = store->elements->datatarray().cbegin(); e != store->elements->datatarray().cend(); ++e) {
		_writer.insert(3, _ecode.data() + 3 * (int)info::mesh->elements->epointers->datatarray()[*e]->code);
	}
	_writer.groupData();

	return store->distribution.process.totalSize;
}

esint VTKLegacy::insertElements(const BoundaryRegionStore *store, const std::vector<char, initless_allocator<char> > &data)
{
	if (store->dimension) {
		esint enodes = 0;
		for (size_t i = 0; i < store->distribution.code.size(); i++) {
			if (store->distribution.code[i].size) {
				enodes += store->distribution.code[i].size * Mesh::edata[i].nodes;
			}
		}

		if (isRoot()) {
			_writer.cells(store->distribution.process.totalSize, store->distribution.process.totalSize + enodes);
		}
		int intsize = 11;
		for (auto e = store->elements->cbegin(); e != store->elements->cend(); ++e) {
			_writer.insert(e->size() > 9 ? 2 : 1, _esize.data() + e->size() * 2 - 2);
			for (auto n = e->begin(); n != e->end(); ++n) {
				esint p = std::lower_bound(store->nodes->datatarray().begin(), store->nodes->datatarray().end(), *n) - store->nodes->datatarray().begin();
				int size = *reinterpret_cast<const int*>(data.data() + (intsize + sizeof(int)) * p + intsize);
				_writer.insert(size, data.data() + (intsize + sizeof(int)) * p);
			}
			_writer.push('\n');
		}
		_writer.groupData();

		if (isRoot()) {
			_writer.celltypes(store->distribution.process.totalSize);
		}
		for (esint e = 0; e < store->distribution.process.size; ++e) {
			_writer.insert(3, _ecode.data() + 3 * (int)store->epointers->datatarray()[e]->code);
		}
		_writer.groupData();

		return store->distribution.process.totalSize;
	} else {
		if (isRoot()) {
			_writer.cells(store->nodeInfo.totalSize, 2 * store->nodeInfo.totalSize);
		}
		for (esint n = 0, i = store->nodeInfo.offset; n < store->nodeInfo.size; ++n, ++i) {
			_writer.cell(1, &i);
		}
		_writer.groupData();

		if (isRoot()) {
			_writer.celltypes(store->nodeInfo.totalSize);
		}
		for (auto n = 0; n < store->nodeInfo.size; ++n) {
			_writer.insert(3, _ecode.data() + 3 * (int)Element::CODE::POINT1);
		}
		_writer.groupData();

		return store->nodeInfo.totalSize;
	}
}

void VTKLegacy::insertData(const step::Step &step, NamedData *data, esint nindices, esint *indices)
{
	if (!storeData(data, step)) {
		return;
	}

	if (data->dataType == NamedData::DataType::SCALAR) {
		for (int d = 0; d < data->dimension; ++d) {
			if (isRoot()) {
				if (data->dimension > 1) {
					_writer.data("SCALARS", data->name + data->coordinateSuffixes[d], "float 1");
				} else {
					_writer.data("SCALARS", data->name, "float 1");
				}
				_writer.description("LOOKUP_TABLE default\n");
			}
			for (esint n = 0; n < nindices; ++n) {
				_writer.float32ln(data->store[indices[n] * data->dimension + d]);
			}
		}
	}

	if (data->dataType == NamedData::DataType::NUMBERED) {
		for (int d = 0; d < data->dimension; ++d) {
			if (isRoot()) {
				_writer.data("SCALARS", data->name + data->numberSuffixes[d], "float 1");
				_writer.description("LOOKUP_TABLE default\n");
			}
			for (esint n = 0; n < nindices; ++n) {
				_writer.float32ln(data->store[indices[n] * data->dimension + d]);
			}
		}
	}

	if (data->dataType == NamedData::DataType::VECTOR) {
		if (isRoot()) {
			_writer.data("VECTORS", data->name, "float");
		}
		for (esint n = 0; n < nindices; ++n) {
			_writer.float32s (                      data->store[indices[n] * data->dimension]);
			_writer.float32s (data->dimension > 1 ? data->store[indices[n] * data->dimension + 1] : .0);
			_writer.float32ln(data->dimension > 2 ? data->store[indices[n] * data->dimension + 2] : .0);
		}
	}

	if (data->dataType == NamedData::DataType::TENSOR_ASYM) {
		if (isRoot()) {
			_writer.data("TENSORS", data->name, "float");
		}
		for (esint n = 0; n < nindices; ++n) {
			_writer.float32s (data->store[indices[n] * data->dimension + 0]);
			_writer.float32s (data->store[indices[n] * data->dimension + 3]);
			_writer.float32s (data->store[indices[n] * data->dimension + 5]);
			_writer.float32s (data->store[indices[n] * data->dimension + 3]);
			_writer.float32s (data->store[indices[n] * data->dimension + 1]);
			_writer.float32s (data->store[indices[n] * data->dimension + 4]);
			_writer.float32s (data->store[indices[n] * data->dimension + 5]);
			_writer.float32s (data->store[indices[n] * data->dimension + 4]);
			_writer.float32ln(data->store[indices[n] * data->dimension + 2]);
		}
	}

	if (data->dataType == NamedData::DataType::TENSOR_SYMM) {
		if (isRoot()) {
			_writer.data("TENSORS", data->name, "float");
		}
		for (esint n = 0; n < nindices; ++n) {
			_writer.float32s (data->store[indices[n] * data->dimension + 0]);
			_writer.float32s (data->store[indices[n] * data->dimension + 3]);
			_writer.float32s (data->store[indices[n] * data->dimension + 5]);
			_writer.float32s (data->store[indices[n] * data->dimension + 6]);
			_writer.float32s (data->store[indices[n] * data->dimension + 1]);
			_writer.float32s (data->store[indices[n] * data->dimension + 4]);
			_writer.float32s (data->store[indices[n] * data->dimension + 8]);
			_writer.float32s (data->store[indices[n] * data->dimension + 7]);
			_writer.float32ln(data->store[indices[n] * data->dimension + 2]);
		}
	}
	_writer.groupData();
}

void VTKLegacy::insertDecomposition(const ElementsRegionStore *store)
{
	auto iterate = [&] (const std::string &name, std::function<double(const ElementsInterval &interval, esint eindex)> callback) {
		if (isRoot()) {
			_writer.data("SCALARS", name, "int 1");
			_writer.description("LOOKUP_TABLE default\n");
		}
		for (auto i = store->eintervals.begin(); i != store->eintervals.end(); ++i) {
			for (esint e = i->begin; e < i->end; ++e) {
				_writer.int32ln(callback(*i, e));
			}
		}
		_writer.groupData();
	};

	iterate("DOMAIN", [&] (const ElementsInterval &interval, esint eindex) {
		return interval.domain;
	});
	esint cluster = info::mesh->clusters->offset;
	Communication::exscan(cluster, MPITools::asynchronous);
	iterate("CLUSTER", [&] (const ElementsInterval &interval, esint eindex) {
		return info::mesh->domains->cluster[interval.domain - info::mesh->domains->offset] + cluster;
	});
	iterate("MPI", [&] (const ElementsInterval &interval, esint eindex) {
		return info::mpi::rank;
	});
}

