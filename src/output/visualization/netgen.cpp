
#include "netgen.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

Netgen::Netgen()
{

}

Netgen::~Netgen()
{

}

void Netgen::updateMesh()
{
	for (size_t i = 0; i < info::mesh->elements->processPerCode.size(); ++i) {
		if (info::mesh->elements->processPerCode[i].totalSize && (int)i != (int)Element::CODE::TETRA4) {
			eslog::error("Netgen writer error: only tetrahedral geometry can be stored in Netgen neutral format.\n");
		}
	}

	if (Visualization::isRoot()) {
		_writer.int32ln(info::mesh->nodes->uniqInfo.totalSize);
	}
	for (esint n = 0, i = info::mesh->nodes->uniqInfo.nhalo; n < info::mesh->nodes->uniqInfo.size; ++n, ++i) {
		_writer.float32(info::mesh->nodes->coordinates->datatarray()[i].x);
		_writer.float32(info::mesh->nodes->coordinates->datatarray()[i].y);
		_writer.float32ln(info::mesh->nodes->coordinates->datatarray()[i].z);
	}
	_writer.groupData();

	esint nelements = 0;
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r) {
		nelements +=  info::mesh->elementsRegions[r]->process.totalSize;
	}
	if (Visualization::isRoot()) {
		_writer.int32ln(nelements);
	}
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r) {
		const ElementsRegionStore *region = info::mesh->elementsRegions[r];

		for (auto e = region->elements->datatarray().begin(); e != region->elements->datatarray().end(); ++e) {
			auto element = (info::mesh->elements->nodes->cbegin() + *e)->data();
			_writer.int32s(r);
			_writer.int32s(info::mesh->nodes->uniqInfo.position[element[0]] + 1);
			_writer.int32s(info::mesh->nodes->uniqInfo.position[element[1]] + 1);
			_writer.int32s(info::mesh->nodes->uniqInfo.position[element[2]] + 1);
			_writer.int32ln(info::mesh->nodes->uniqInfo.position[element[3]] + 1);
		}
	}
	_writer.groupData();

	esint nboundary = 0;
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			nboundary +=  info::mesh->boundaryRegions[r]->process.totalSize;
		}
	}
	if (Visualization::isRoot()) {
		_writer.int32ln(nboundary);
	}
	for (size_t r = 1, bindex = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];

		if (region->dimension) {
			for (auto e = region->procNodes->begin(); e != region->procNodes->end(); ++e) {
				auto element = e->data();
				_writer.int32s(bindex);
				_writer.int32s(info::mesh->nodes->uniqInfo.position[element[0]] + 1);
				_writer.int32s(info::mesh->nodes->uniqInfo.position[element[1]] + 1);
				_writer.int32ln(info::mesh->nodes->uniqInfo.position[element[2]] + 1);
			}
			++bindex;
		}
	}
	_writer.groupData();

	_writer.commitFile(_path + _name + ".mesh");
	_writer.reorder();
	_writer.write();
}


void Netgen::updateSolution()
{
	// TODO
}
