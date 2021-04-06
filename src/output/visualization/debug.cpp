
#include "debug.h"
#include "visualization.h"
#include "writer/vtkwritter.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/communication.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "mesh/mesh.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/contactstore.h"

#include <cstddef>
#include <fstream>
#include <vector>

using namespace espreso;

DebugOutput::DebugOutput(double clusterShrinkRatio, double domainShrinkRatio, bool withDomains)
: _mesh(*info::mesh), _ccenter(Point()), _dcenters(NULL),
  _clusterShrinkRatio(clusterShrinkRatio), _domainShrinkRatio(domainShrinkRatio)
{
	_path = utils::createDirectory({ std::string(eslog::path()), "DEBUG_VISUALIZATION" });

	if (Visualization::isRoot()) {
		_writer.description("# vtk DataFile Version 2.0\n");
		_writer.description("EXAMPLE\n");
		_writer.description("ASCII\n");
		_writer.description("DATASET UNSTRUCTURED_GRID\n");
	}

	if (withDomains) {
		if (_mesh.nodes->domains == NULL) {
			mesh::computeNodeDomainDistribution();
		}

		std::vector<Point> dcenters(_mesh.elements->ndomains);
		std::vector<esint> dcounter(_mesh.elements->ndomains);
		auto domains = _mesh.nodes->domains->begin();
		for (esint n = 0; n < _mesh.nodes->size; ++n, ++domains) {
			for (auto d = domains->begin(); d != domains->end(); ++d) {
				if (_mesh.elements->firstDomain <= *d && *d < _mesh.elements->firstDomain + _mesh.elements->ndomains) {
					dcenters[*d - _mesh.elements->firstDomain] += _mesh.nodes->coordinates->datatarray()[n];
					dcounter[*d - _mesh.elements->firstDomain] += 1;
				}
			}
		}

		_dcenters = new Point[_mesh.elements->ndomains];
		for (size_t d = 0; d < dcenters.size(); ++d) {
			_dcenters[d] = dcenters[d] / dcounter[d];
		}
	}

	for (esint n = 0; n < _mesh.nodes->size; ++n) {
		_ccenter += _mesh.nodes->coordinates->datatarray()[n];
	}
	_ccenter /= _mesh.nodes->size;
}

void DebugOutput::points(esint nother, esint &noffset, esint &nsize)
{
	noffset = _mesh.nodes->size;
	nsize = Communication::exscan(noffset);

	if (Visualization::isRoot()) {
		_writer.points(nsize + nother);
	}
	for (esint n = 0; n < _mesh.nodes->size; ++n) {
		Point p = Visualization::shrink(_mesh.nodes->coordinates->datatarray()[n], _ccenter, Point(), _clusterShrinkRatio, 1);
		_writer.point(p.x, p.y, p.z);
	}
	_writer.groupData();
}

void DebugOutput::pointsInDomains(esint nother, esint &noffset, esint &nsize)
{
	noffset = _mesh.nodes->domains->datatarray().size();
	nsize = Communication::exscan(noffset);

	esint gnother;
	Communication::allReduce(&nother, &gnother, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		_writer.points(nsize + gnother);
	}
	auto c = _mesh.nodes->coordinates->datatarray().begin();
	for (auto n = _mesh.nodes->domains->begin(); n != _mesh.nodes->domains->end(); ++n, ++c) {
		for (auto d = n->begin(); d != n->end(); ++d) {
			Point p = Visualization::shrink(*c, _ccenter, _dcenters[*d - _mesh.elements->firstDomain], _clusterShrinkRatio, _domainShrinkRatio);
			_writer.point(p.x, p.y, p.z);
		}
	}
	_writer.groupData();
}

esint DebugOutput::elements(esint noffset, esint nother, esint nothernodes)
{
	esint esize = _mesh.elementsRegions.front()->elements->structures(), gesize;
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	esint ensize = _mesh.elements->procNodes->datatarray().size(), gensize;
	Communication::allReduce(&ensize, &gensize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		_writer.cells(gesize + nother, gesize + gensize + nother + nothernodes);
	}

	for (auto e = _mesh.elements->procNodes->begin(); e != _mesh.elements->procNodes->end(); ++e) {
		_writer.cell(e->size(), e->data(), noffset);
	}
	_writer.groupData();

	return gesize;
}

esint DebugOutput::elementsInDomains(esint noffset, esint nother, esint nothernodes)
{
	esint esize = _mesh.elementsRegions.front()->elements->structures(), gesize;
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	esint ensize = _mesh.elements->procNodes->datatarray().size(), gensize;
	Communication::allReduce(&ensize, &gensize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		_writer.cells(gesize, gesize + gensize);
	}

	for (auto e = _mesh.elements->domainNodes->begin(); e != _mesh.elements->domainNodes->end(); ++e) {
		_writer.cell(e->size(), e->data(), noffset);
	}
	_writer.groupData();

	return esize;
}

void DebugOutput::etypes(esint esize, esint nother)
{
	if (Visualization::isRoot()) {
		_writer.celltypes(esize + nother);
	}
	for (auto e = _mesh.elements->epointers->datatarray().begin(); e != _mesh.elements->epointers->datatarray().end(); ++e) {
		_writer.type((*e)->code);
	}
	_writer.groupData();
}

DebugOutput::~DebugOutput()
{
	if (_dcenters) { delete[] _dcenters; }
}

void DebugOutput::mesh(double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}

	DebugOutput output(clusterShrinkRatio, 1, false);

	esint noffset, nsize;
	output.points(0, noffset, nsize);
	output._writer.groupData();

	esint esize = output.elements(noffset, 0, 0);
	output._writer.groupData();

	output.etypes(esize, 0);
	output._writer.groupData();

	output._writer.commitFile(output._path + "mesh.vtk");
	output._writer.reorder();
	output._writer.write();
}

void DebugOutput::faceNeighbors()
{
	if (!info::ecf->output.debug) {
		return;
	}

	DebugOutput output(1, 1, false);

	if (output._mesh.elements->centers == NULL) {
		mesh::computeElementsCenters();
	}

	esint psize = output._mesh.elements->centers->datatarray().size(), gpsize;
	Communication::allReduce(&psize, &gpsize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	if (Visualization::isRoot()) {
		output._writer.points(gpsize);
	}
	for (size_t n = 0; n < output._mesh.elements->centers->datatarray().size(); ++n) {
		Point &p = output._mesh.elements->centers->datatarray()[n];
		output._writer.point(p.x, p.y, p.z);
	}
	output._writer.groupData();

	esint esize = 0, gesize;
	auto id = output._mesh.elements->IDs->begin();
	for (auto neigh = output._mesh.elements->faceNeighbors->begin(); neigh != output._mesh.elements->faceNeighbors->end(); ++neigh, ++id) {
		for (auto n = neigh->begin(); n != neigh->end(); ++n) {
			if (id->front() < *n) {
				++esize;
			}
		}
	}
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		output._writer.cells(gesize, gesize + 2 * gesize);
	}

	id = output._mesh.elements->IDs->begin();
	for (auto neigh = output._mesh.elements->faceNeighbors->begin(); neigh != output._mesh.elements->faceNeighbors->end(); ++neigh, ++id) {
		for (auto n = neigh->begin(); n != neigh->end(); ++n) {
			if (id->front() < *n) {
				esint nn[2] = { id->front(), *n };
				output._writer.cell(2, nn);
			}
		}
	}
	output._writer.groupData();

	if (Visualization::isRoot()) {
		output._writer.celltypes(gesize);
	}
	id = output._mesh.elements->IDs->begin();
	for (auto neigh = output._mesh.elements->faceNeighbors->begin(); neigh != output._mesh.elements->faceNeighbors->end(); ++neigh, ++id) {
		for (auto n = neigh->begin(); n != neigh->end(); ++n) {
			if (id->front() < *n) {
				output._writer.type(Element::CODE::LINE2);
			}
		}
	}
	output._writer.groupData();

	output._writer.commitFile(output._path + "faceNeighbors.vtk");
	output._writer.reorder();
	output._writer.write();
}

void DebugOutput::meshDual(std::vector<esint> &frames, std::vector<esint> &neighbors)
{
	if (!info::ecf->output.debug) {
		return;
	}

	DebugOutput output(1, 1, false);

	if (output._mesh.elements->centers == NULL) {
		mesh::computeElementsCenters();
	}

	esint psize = output._mesh.elements->centers->datatarray().size(), gpsize;
	Communication::allReduce(&psize, &gpsize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	if (Visualization::isRoot()) {
		output._writer.points(gpsize);
	}
	for (size_t n = 0; n < output._mesh.elements->centers->datatarray().size(); ++n) {
		Point &p = output._mesh.elements->centers->datatarray()[n];
		output._writer.point(p.x, p.y, p.z);
	}
	output._writer.groupData();

	esint esize = 0, gesize;
	auto id = output._mesh.elements->IDs->begin();
	for (size_t e = 0; e < frames.size() - 1; ++e, ++id) {
		for (auto n = frames[e]; n != frames[e + 1]; ++n) {
			if (id->front() < neighbors[n]) {
				++esize;
			}
		}
	}
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		output._writer.cells(gesize, gesize + 2 * gesize);
	}

	id = output._mesh.elements->IDs->begin();
	for (size_t e = 0; e < frames.size() - 1; ++e, ++id) {
		for (auto n = frames[e]; n != frames[e + 1]; ++n) {
			if (id->front() < neighbors[n]) {
				esint nn[2] = { id->front(), neighbors[n] };
				output._writer.cell(2, nn);
			}
		}
	}
	output._writer.groupData();

	if (Visualization::isRoot()) {
		output._writer.celltypes(gesize);
	}
	id = output._mesh.elements->IDs->begin();
	for (size_t e = 0; e < frames.size() - 1; ++e, ++id) {
		for (auto n = frames[e]; n != frames[e + 1]; ++n) {
			if (id->front() < neighbors[n]) {
				output._writer.type(Element::CODE::LINE2);
			}
		}
	}
	output._writer.groupData();

	output._writer.commitFile(output._path + "meshDual.vtk");
	output._writer.reorder();
	output._writer.write();
}

void DebugOutput::corners(double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}
}

void DebugOutput::innerFixPoints(double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}
}

void DebugOutput::surfaceFixPoints(double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}
}

void DebugOutput::contact(double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}

	DebugOutput output(clusterShrinkRatio, 1, false);

	esint cells = 0, gcells, points = 0, poffset, noffset, nsize;;
	if (output._mesh.contacts->intersections != NULL) {
		cells += output._mesh.contacts->intersections->datatarray().size();
		points += 3 * output._mesh.contacts->intersections->datatarray().size();
	}
	poffset = points;
	points = Communication::exscan(poffset);
	Communication::allReduce(&cells, &gcells, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	output.points(points, noffset, nsize);
	if (output._mesh.contacts->intersections != NULL) {
		for (size_t i = 0; i < output._mesh.contacts->intersections->datatarray().size(); ++i) {
			for (int p = 0; p < 3; ++p) {
				Point pp = Visualization::shrink(output._mesh.contacts->intersections->datatarray()[i].p[p], output._ccenter, Point(), .95, 1);
				output._writer.point(pp.x, pp.y, pp.z);
			}
		}
	}
	output._writer.groupData();

	esint esize = output.elements(noffset, gcells, points);
	if (output._mesh.contacts->intersections != NULL) {
		poffset += nsize;
		for (size_t i = 0; i < output._mesh.contacts->intersections->datatarray().size(); ++i) {
			esint nn[3] = { poffset++, poffset++, poffset++};
			output._writer.cell(3, nn);
		}
	}
	output._writer.groupData();

	if (output._mesh.contacts->intersections != NULL) {
		output.etypes(esize, gcells);
		for (size_t i = 0; i < output._mesh.contacts->intersections->datatarray().size(); ++i) {
			output._writer.type(Element::CODE::TRIANGLE3);
		}
	}
	output._writer.groupData();

	if (Visualization::isRoot()) {
		output._writer.celldata(esize + gcells);
		output._writer.data("SCALARS", "MESH", "int 1");
		output._writer.description("LOOKUP_TABLE default\n");
	}
	for (esint e = 0; e < output._mesh.elements->size; ++e) {
		output._writer.int32ln(1);
	}
	output._writer.groupData();

	for (esint e = 0; e < cells; ++e) {
		output._writer.int32ln(0);
	}
	output._writer.groupData();

	output._writer.commitFile(output._path + "contact.vtk");
	output._writer.reorder();
	output._writer.write();
}

void DebugOutput::surface(const char* name, double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}

	DebugOutput output(clusterShrinkRatio, domainShrinkRatio, false);

	esint points = 0, noffset, nsize;
	output.points(points, noffset, nsize);

	// elements

	esint gesize, esize = output._mesh.surface->epointers->datatarray().size();
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	esint gensize, ensize = output._mesh.surface->enodes->datatarray().size();
	Communication::allReduce(&ensize, &gensize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		output._writer.cells(gesize, gesize + gensize);
	}

	for (auto e = output._mesh.surface->enodes->begin(); e != output._mesh.surface->enodes->end(); ++e) {
		output._writer.cell(e->size(), e->data(), noffset);
	}
	output._writer.groupData();

	if (Visualization::isRoot()) {
		output._writer.celltypes(gesize);
	}
	for (auto e = output._mesh.surface->epointers->datatarray().begin(); e != output._mesh.surface->epointers->datatarray().end(); ++e) {
		output._writer.type((*e)->code);
	}
	output._writer.groupData();
	output._writer.commitFile(output._path + std::string(name) + ".vtk");
	output._writer.reorder();
	output._writer.write();
}

void DebugOutput::closeElements(double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}
	if (info::mesh->contacts == NULL || info::mesh->contacts->closeElements == NULL) {
		return;
	}

	Point center;
	for (auto n = info::mesh->nodes->coordinates->datatarray().begin(); n != info::mesh->nodes->coordinates->datatarray().end(); ++n) {
		center += *n;
	}
	center /= info::mesh->nodes->size;

	std::string path = utils::createDirectory({ std::string(eslog::path()), "DEBUG_VISUALIZATION" });
	std::ofstream os(path + "/closeElements" + std::to_string(info::mpi::rank) + ".vtk");

	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	size_t points = info::mesh->contacts->elements->structures() + info::mesh->contacts->closeElements->datatarray().size();
	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		points += info::mesh->contacts->gnsurface[n].size() + info::mesh->contacts->gncloseElements[n]->datatarray().size();
	}

	std::vector<Point> centers(info::mpi::size);
	Communication::allGather(&center, centers.data(), 3, MPI_DOUBLE);

	os << "POINTS " << points << " float\n";
	auto closest = info::mesh->contacts->closeElements->cbegin();
	for (auto e = info::mesh->contacts->elements->cbegin(); e !=info::mesh->contacts->elements->cend(); ++e, ++closest) {
		Point center;
		for (auto n = e->begin(); n != e->end(); ++n) {
			center += *n;
		}
		center /= e->size();
		center = Visualization::shrink(center, center, center, clusterShrinkRatio, 1);
		os << center.x << " " << center.y << " " << center.z << "\n";

		for (auto ne = closest->begin(); ne != closest->end(); ++ne) {
			auto ce = info::mesh->contacts->elements->cbegin() + *ne;

			Point ncenter;
			for (auto n = ce->begin(); n != ce->end(); ++n) {
				ncenter += *n;
			}
			ncenter /= ce->size();

			ncenter = Visualization::shrink(ncenter, center, ncenter, clusterShrinkRatio, 1);
			ncenter = center + (ncenter - center) / 2;
			os << ncenter.x << " " << ncenter.y << " " << ncenter.z << "\n";
		}
	}

	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		auto closest = info::mesh->contacts->gncloseElements[n]->cbegin();
		for (size_t i = 0; i < info::mesh->contacts->gnsurface[n].size(); ++i, ++closest) {
			auto e = info::mesh->contacts->elements->cbegin() + info::mesh->contacts->gnsurface[n][i];
			Point center;
			for (auto nn = e->begin(); nn != e->end(); ++nn) {
				center += *nn;
			}
			center /= e->size();
			center = Visualization::shrink(center, center, center, clusterShrinkRatio, 1);
			os << center.x << " " << center.y << " " << center.z << "\n";

			for (auto ne = closest->begin(); ne != closest->end(); ++ne) {
				auto ce = info::mesh->contacts->gnecoords[n]->cbegin() + *ne;

				Point ncenter;
				for (auto nn = ce->begin(); nn != ce->end(); ++nn) {
					ncenter += *nn;
				}
				ncenter /= ce->size();

				ncenter = Visualization::shrink(ncenter, centers[info::mesh->contacts->gneighbors[n]], ncenter, clusterShrinkRatio, 1);
				ncenter = center + (ncenter - center) / 2;
				os << ncenter.x << " " << ncenter.y << " " << ncenter.z << "\n";
			}
		}
	}

	os << "\n";

	size_t cells = info::mesh->contacts->closeElements->datatarray().size();
	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		cells += info::mesh->contacts->gncloseElements[n]->datatarray().size();
	}

	os << "CELLS " << cells << " " << 3 * cells << "\n";
	esint eindex = 0, noffset;
	for (auto e = info::mesh->contacts->closeElements->cbegin(); e != info::mesh->contacts->closeElements->cend(); ++e) {
		noffset = 1;
		for (auto n = e->begin(); n != e->end(); ++n, ++noffset) {
			os << "2 " << eindex << " " << eindex + noffset << "\n";
		}
		eindex += noffset;
	}
	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		for (auto e = info::mesh->contacts->gncloseElements[n]->cbegin(); e != info::mesh->contacts->gncloseElements[n]->cend(); ++e) {
			noffset = 1;
			for (auto n = e->begin(); n != e->end(); ++n, ++noffset) {
				os << "2 " << eindex << " " << eindex + noffset << "\n";
			}
			eindex += noffset;
		}
	}
	os << "\n";

	os << "CELL_TYPES " << cells << "\n";
	Element::CODE ecode = Element::CODE::LINE2;
	for (size_t n = 0; n < cells; ++n) {
		os << VTKASCIIWritter::ecode(ecode) << "\n";
	}
	os << "\n";
}
