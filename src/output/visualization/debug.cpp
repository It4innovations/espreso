
#include "debug.h"
#include "visualization.h"
#include "writer/vtkwritter.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/sysutils.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "mesh/mesh.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
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
	_path = utils::createDirectory({ info::ecf->outpath, "DEBUG_VISUALIZATION" });

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

		std::vector<Point> dcenters(_mesh.domains->size);
		std::vector<esint> dcounter(_mesh.domains->size);
		auto domains = _mesh.nodes->domains->begin();
		for (esint n = 0; n < _mesh.nodes->size; ++n, ++domains) {
			for (auto d = domains->begin(); d != domains->end(); ++d) {
				if (_mesh.domains->offset <= *d && *d < _mesh.domains->offset + _mesh.domains->size) {
					dcenters[*d - _mesh.domains->offset] += _mesh.nodes->coordinates->datatarray()[n];
					dcounter[*d - _mesh.domains->offset] += 1;
				}
			}
		}

		_dcenters = new Point[_mesh.domains->size];
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
			Point p = Visualization::shrink(*c, _ccenter, _dcenters[*d - _mesh.domains->offset], _clusterShrinkRatio, _domainShrinkRatio);
			_writer.point(p.x, p.y, p.z);
		}
	}
	_writer.groupData();
}

esint DebugOutput::elements(esint noffset, esint nother, esint nothernodes)
{
	esint esize = _mesh.elementsRegions.front()->elements->structures(), gesize;
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	esint ensize = _mesh.elements->nodes->datatarray().size(), gensize;
	Communication::allReduce(&ensize, &gensize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		_writer.cells(gesize + nother, gesize + gensize + nother + nothernodes);
	}

	for (auto e = _mesh.elements->nodes->begin(); e != _mesh.elements->nodes->end(); ++e) {
		_writer.cell(e->size(), e->data(), noffset);
	}
	_writer.groupData();

	return gesize;
}

esint DebugOutput::elementsInDomains(esint noffset, esint nother, esint nothernodes)
{
	esint esize = _mesh.elementsRegions.front()->elements->structures(), gesize;
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	esint ensize = _mesh.elements->nodes->datatarray().size(), gensize;
	Communication::allReduce(&ensize, &gensize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		_writer.cells(gesize, gesize + gensize);
	}

	for (auto e = _mesh.domains->nodes->begin(); e != _mesh.domains->nodes->end(); ++e) {
		_writer.cell(e->size(), e->data(), noffset);
	}
	_writer.groupData();

	return esize;
}

void DebugOutput::data(const std::string &name, const std::vector<Point> &points, const std::vector<std::vector<esint> > &cells, const std::vector<esint> &celltypes, const std::vector<std::vector<double> > &celldata)
{
	std::ofstream os(name + ".vtk");
	os << "# vtk DataFile Version 2.0\n";
	os << "CLIP\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	os << "POINTS " << points.size() << " float\n";
	for (size_t pp = 0; pp < points.size(); ++pp) {
		os << points[pp].x << " " << points[pp].y << " " << points[pp].z << "\n";
	}
	os << "\n";

	esint cnodes = 0;
	for (size_t cc = 0; cc < cells.size(); ++cc) {
		cnodes += cells[cc].size();
	}

	os << "CELLS " << cells.size() << " " << cells.size() + cnodes << "\n";
	for (size_t cc = 0; cc < cells.size(); ++cc) {
		os << cells[cc].size();
		for (size_t cn = 0; cn < cells[cc].size(); ++cn) {
			os << " " << cells[cc][cn];
		}
		os << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << celltypes.size() << "\n";
	for (size_t cc = 0; cc < celltypes.size(); ++cc) {
		os << celltypes[cc] << "\n";
	}
	os << "\n";

	os << "CELL_DATA " << cells.size() << "\n";
	os << "SCALARS DISTANCE float 1\n";
	os << "LOOKUP_TABLE default\n";
	for (size_t dd = 0; dd < celldata.size(); ++dd) {
		for (size_t cd = 0; cd < celldata[dd].size(); ++cd) {
			os << celldata[dd][cd] << "\n";
		}
	}
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
	if (!info::ecf->output.debug || info::mesh->contacts->sparseSide == NULL) {
		return;
	}

	DebugOutput output(clusterShrinkRatio, domainShrinkRatio, false);

	auto *sside = info::mesh->contacts->sparseSide;
	auto *dside = info::mesh->contacts->denseSide;

	esint cells = 0, gcells, points = 0, poffset;
	for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
		for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
			if (!d->skip) {
				points += 3 * d->triangles;
				cells += d->triangles;
			}
		}
	}

	poffset = points;
	points = Communication::exscan(poffset);
	Communication::allReduce(&cells, &gcells, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		output._writer.points(points);
	}

	esint triangle = 0;
	for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
		for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
			if (!d->skip) {
				for (esint t = 0; t < d->triangles; ++t) {
					const Triangle &tr = output._mesh.contacts->intersections->datatarray()[triangle++];
					output._writer.point(tr.p[0].x, tr.p[0].y, tr.p[0].z);
					output._writer.point(tr.p[1].x, tr.p[1].y, tr.p[1].z);
					output._writer.point(tr.p[2].x, tr.p[2].y, tr.p[2].z);
				}
			} else {
				triangle += d->triangles;
			}
		}
	}
	output._writer.groupData();

	if (Visualization::isRoot()) {
		output._writer.cells(gcells, 4 * gcells);
	}
	for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
		for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
			if (!d->skip) {
				for (esint t = d->triangleOffset; t < d->triangleOffset + d->triangles; ++t) {
					esint nn[3] = { poffset++, poffset++, poffset++};
					output._writer.cell(3, nn);
				}
			}
		}
	}
	output._writer.groupData();

	if (Visualization::isRoot()) {
		output._writer.celltypes(gcells);
	}
	for (esint i = 0; i < cells; ++i) {
		output._writer.type(Element::CODE::TRIANGLE3);
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
		output._writer.cell(e->size(), e->data(), output._mesh.surface->nodes->datatarray().data(), noffset);
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

void DebugOutput::warpedNormals(const char* name, double clusterShrinkRatio, double domainShrinkRatio)
{
	if (!info::ecf->output.debug) {
		return;
	}

	DebugOutput output(1, 1, false);

	SurfaceStore *surf = output._mesh.surface;

	esint psize = surf->enodes->datatarray().size() + 2 * surf->size, gpsize;
	Communication::allReduce(&psize, &gpsize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);
	if (Visualization::isRoot()) {
		output._writer.points(gpsize);
	}
	auto enodes = surf->enodes->begin();
	auto normal = surf->normal->datatarray().begin();
	auto center = surf->base->datatarray().begin();
	for (esint e = 0; e < surf->size; ++e, ++enodes, ++normal, ++center) {
		Point &p = *center;
		Point &n = *normal;
		double scale = 0.001;
		output._writer.point(p.x, p.y, p.z);
		output._writer.point(p.x + scale * n.x, p.y + scale * n.y, p.z + scale * n.z);
		for (auto nn = enodes->begin(); nn != enodes->end(); ++nn) {
			Point pn = surf->coordinates->datatarray()[*nn];
			double distance = (pn - p) * n;
			pn = pn - n * distance;
			output._writer.point(pn.x, pn.y, pn.z);
		}
	}
	output._writer.groupData();

	esint esize = 2 * surf->size, gesize;
	Communication::allReduce(&esize, &gesize, 1, MPITools::getType<esint>().mpitype, MPI_SUM);

	if (Visualization::isRoot()) {
		output._writer.cells(gesize, gesize + gpsize);
	}
	enodes = surf->enodes->begin();
	for (esint e = 0, n = 0; e < surf->size; ++e, ++enodes) {
		output._writer.int32s(2);
		output._writer.int32s(n++);
		output._writer.int32ln(n++);

		output._writer.int32s(enodes->size());
		for (auto nn = enodes->begin(); nn != enodes->end(); ++nn) {
			output._writer.int32s(n++);
		}
		output._writer.push('\n');
	}

	output._writer.groupData();

	if (Visualization::isRoot()) {
		output._writer.celltypes(gesize);
	}
	auto epointers = surf->epointers->datatarray().begin();
	for (esint e = 0; e < surf->size; ++e, ++epointers) {
		output._writer.type(Element::CODE::LINE2);
		output._writer.type((*epointers)->code);
	}
	output._writer.groupData();

	output._writer.commitFile(output._path + std::string(name) + ".vtk");
	output._writer.reorder();
	output._writer.write();
}

void DebugOutput::closeElements(double clusterShrinkRatio, double domainShrinkRatio)
{
//	if (!info::ecf->output.debug) {
//		return;
//	}
//	if (info::mesh->contacts == NULL || info::mesh->contacts->localPairs == NULL) {
//		return;
//	}
//
//	Point center;
//	for (auto n = info::mesh->nodes->coordinates->datatarray().begin(); n != info::mesh->nodes->coordinates->datatarray().end(); ++n) {
//		center += *n;
//	}
//	center /= info::mesh->nodes->size;
//
//	std::string path = utils::createDirectory({ std::string(eslog::path()), "DEBUG_VISUALIZATION" });
//	std::ofstream os(path + "/closeElements" + std::to_string(info::mpi::rank) + ".vtk");
//
//	os << "# vtk DataFile Version 2.0\n";
//	os << "EXAMPLE\n";
//	os << "ASCII\n";
//	os << "DATASET UNSTRUCTURED_GRID\n\n";
//
//	size_t points = info::mesh->surface->enodes->structures() + info::mesh->contacts->localPairs->datatarray().size() + info::mesh->contacts->neighPairs->datatarray().size() / 2;
//
//	os << "POINTS " << points << " float\n";
//	std::vector<float> distance;
//	auto enodes = info::mesh->surface->enodes->cbegin();
//	auto local = info::mesh->contacts->localPairs->cbegin();
//	auto neigh = info::mesh->contacts->neighPairs->cbegin();
//	for (esint e = 0; e < info::mesh->surface->size; ++e, ++enodes, ++local, ++neigh) {
//		Point center;
//		for (auto n = enodes->begin(); n != enodes->end(); ++n) {
//			center += info::mesh->surface->coordinates->datatarray()[*n];
//		}
//		center /= enodes->size();
//		os << center.x << " " << center.y << " " << center.z << "\n";
//		for (auto ee = local->begin(); ee != local->end(); ++ee) {
//			Point tocenter;
//			auto other = info::mesh->surface->enodes->cbegin() + *ee;
//			for (auto n = other->begin(); n != other->end(); ++n) {
//				tocenter += info::mesh->surface->coordinates->datatarray()[*n];
//			}
//			tocenter /= other->size();
//			os << tocenter.x << " " << tocenter.y << " " << tocenter.z << "\n";
//			distance.push_back((tocenter - center).length());
//		}
//		for (auto ee = neigh->begin(); ee != neigh->end(); ++ee) {
//			Point tocenter;
//			esint nn = *ee++;
//			auto other = info::mesh->contacts->surfaces[nn]->enodes->cbegin() + *ee;
//			for (auto n = other->begin(); n != other->end(); ++n) {
//				tocenter += info::mesh->contacts->surfaces[nn]->coordinates->datatarray()[*n];
//			}
//			tocenter /= other->size();
//			os << tocenter.x << " " << tocenter.y << " " << tocenter.z << "\n";
//			distance.push_back((tocenter - center).length());
//		}
//	}
//	os << "\n";
//
//	size_t cells = info::mesh->contacts->localPairs->datatarray().size() + info::mesh->contacts->neighPairs->datatarray().size() / 2;
//
//	os << "CELLS " << cells << " " << 3 * cells << "\n";
//	esint eindex, noffset = 0;
//	local = info::mesh->contacts->localPairs->cbegin();
//	neigh = info::mesh->contacts->neighPairs->cbegin();
//	for (esint e = 0; e < info::mesh->surface->size; ++e, ++local, ++neigh) {
//		eindex = noffset++;
//		for (auto ee = local->begin(); ee != local->end(); ++ee, ++noffset) {
//			os << "2 " << eindex << " " << noffset << "\n";
//		}
//		for (auto ee = neigh->begin(); ee != neigh->end(); ee += 2, ++noffset) {
//			os << "2 " << eindex << " " << noffset << "\n";
//		}
//	}
//
//	os << "\n";
//
//	os << "CELL_TYPES " << cells << "\n";
//	Element::CODE ecode = Element::CODE::LINE2;
//	for (size_t n = 0; n < cells; ++n) {
//		os << VTKASCIIWritter::ecode(ecode) << "\n";
//	}
//	os << "\n";
//
//	os << "CELL_DATA " << cells << "\n";
//	os << "SCALARS DISTANCE float 1\n";
//	os << "LOOKUP_TABLE default\n";
//	for (size_t n = 0; n < cells; ++n) {
//		os << distance[n] << "\n";
//	}
//	os << "\n";
}
