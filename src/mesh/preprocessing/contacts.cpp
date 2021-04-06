
#include "meshpreprocessing.h"
#include "weileratherton.h"

#include "mesh/element.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/contactstore.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "basis/logging/timelogger.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "output/visualization/debug.h"

#include <algorithm>
#include <numeric>
#include <iomanip>
#include <cfloat>
#include <fstream>

#include "basis/utilities/print.h"

#define MIN_SLAVE_COVER_RATIO 0.3
#define TOLERATED_SLAVE_COVER_RATIO_FOR_DUAL_SHAPE_COEFFICIENTS_ON_WHOLE_ELEMENT 0.3
#define BE_VALUE_TRESHOLD 1e-15

namespace espreso {
namespace mesh {

void computeBodiesSurface()
{
	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointers(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = info::mesh->elements->procNodes->cbegin(t);
		auto neighs = info::mesh->elements->faceNeighbors->cbegin(t);
		const auto &epointers = info::mesh->elements->epointers->datatarray();

		std::vector<esint> fdist, fdata, ecounter((int)Element::CODE::SIZE);
		std::vector<Element*> fpointer;
		if (t == 0) {
			fdist.push_back(0);
		}

		for (size_t e = info::mesh->elements->distribution[t]; e < info::mesh->elements->distribution[t + 1]; ++e, ++neighs, ++nodes) {
			for (size_t n = 0; n < neighs->size(); ++n) {
				if (neighs->at(n) == -1) {
					auto face = epointers[e]->faces->begin() + n;
					for (auto f = face->begin(); f != face->end(); ++f) {
						fdata.push_back(nodes->at(*f));
					}
					fdist.push_back(fdata.size());
					fpointer.push_back(epointers[e]->facepointers->datatarray()[n]);
					++ecounter[(int)fpointer.back()->code];
				}
			}
		}

		facesDistribution[t].swap(fdist);
		faces[t].swap(fdata);
		fpointers[t].swap(fpointer);
		ecounters[t].swap(ecounter);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t e = 0; e < ecounters[0].size(); e++) {
			ecounters[0][e] += ecounters[t][e];
		}
	}

	serializededata<esint, Element*>::balance(1, fpointers);
	info::mesh->surface->epointers = new serializededata<esint, Element*>(1, fpointers);
	info::mesh->surface->ecounters = ecounters[0];

	info::mesh->surface->edistribution = info::mesh->surface->epointers->datatarray().distribution();

	if (info::mesh->surface->ecounters[(int)Element::CODE::TRIANGLE3] == (esint)info::mesh->surface->edistribution.back()) {
		serializededata<esint, esint>::balance(3, faces, &info::mesh->surface->edistribution);
		info::mesh->surface->enodes = new serializededata<esint, esint>(3, faces);
		info::mesh->surface->triangles = info::mesh->surface->enodes;
		info::mesh->surface->tdistribution = info::mesh->surface->edistribution;
	} else {
		utils::threadDistributionToFullDistribution(facesDistribution);
		serializededata<esint, esint>::balance(facesDistribution, faces, &info::mesh->surface->edistribution);
		info::mesh->surface->enodes = new serializededata<esint, esint>(facesDistribution, faces);
	}

	DebugOutput::surface("surface.bodies");

	eslog::checkpointln("MESH: BODY SURFACE COMPUTED");
}

void computeSurfaceLocations()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<Point> > points(threads);
	std::vector<std::vector<esint> > nids(threads);
	std::vector<std::vector<Element*> > epointers(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		points[t].resize(info::mesh->surface->enodes->cbegin(t + 1)->begin() - info::mesh->surface->enodes->cbegin(t)->begin());
		nids[t].resize(info::mesh->surface->enodes->cbegin(t + 1)->begin() - info::mesh->surface->enodes->cbegin(t)->begin());
		epointers[t].insert(epointers[t].end(), info::mesh->surface->epointers->datatarray().begin(t), info::mesh->surface->epointers->datatarray().end(t));
		size_t i = 0;
		for (auto e = info::mesh->surface->enodes->cbegin(t); e != info::mesh->surface->enodes->cend(t); ++e) {
			for (auto n = e->begin(); n != e->end(); ++n, ++i) {
				points[t][i] = info::mesh->nodes->coordinates->datatarray()[*n];
				nids[t][i] = info::mesh->nodes->IDs->datatarray()[*n];
			}
		}
	}

	if (info::mesh->surface->enodes->boundarytarray().size()) {
		info::mesh->contacts->elements = new serializededata<esint, Point>(tarray<esint>(info::mesh->surface->enodes->boundarytarray()), points);
		info::mesh->contacts->enodes = new serializededata<esint, esint>(tarray<esint>(info::mesh->surface->enodes->boundarytarray()), nids);
	} else {
		info::mesh->contacts->elements = new serializededata<esint, Point>(3, points);
		info::mesh->contacts->enodes = new serializededata<esint, esint>(3, nids);
	}
	info::mesh->contacts->epointers = new serializededata<esint, Element*>(1, epointers);

	size_t precision = 0;
	double epsilon = 1e-6;
	while (true) {
		double value = info::mesh->contacts->eps * pow(10, precision);
		if (std::round(value) <= value + epsilon && value - epsilon <= std::round(value)) {
			break;
		} else {
			++precision;
		}
	}

	// 1. compute bounding box
	//////////////////////////

	info::mesh->contacts->boundingBox[0] = info::mesh->nodes->coordinates->datatarray().front();
	info::mesh->contacts->boundingBox[1] = info::mesh->nodes->coordinates->datatarray().front();
	for (esint n = 0; n < info::mesh->nodes->size; ++n) {
		info::mesh->contacts->boundingBox[0].x = std::min(info::mesh->contacts->boundingBox[0].x, info::mesh->nodes->coordinates->datatarray()[n].x);
		info::mesh->contacts->boundingBox[0].y = std::min(info::mesh->contacts->boundingBox[0].y, info::mesh->nodes->coordinates->datatarray()[n].y);
		info::mesh->contacts->boundingBox[0].z = std::min(info::mesh->contacts->boundingBox[0].z, info::mesh->nodes->coordinates->datatarray()[n].z);
		info::mesh->contacts->boundingBox[1].x = std::max(info::mesh->contacts->boundingBox[1].x, info::mesh->nodes->coordinates->datatarray()[n].x);
		info::mesh->contacts->boundingBox[1].y = std::max(info::mesh->contacts->boundingBox[1].y, info::mesh->nodes->coordinates->datatarray()[n].y);
		info::mesh->contacts->boundingBox[1].z = std::max(info::mesh->contacts->boundingBox[1].z, info::mesh->nodes->coordinates->datatarray()[n].z);
	}

	auto rounddown = [&] (double &value) {
		int rounder = info::mesh->contacts->eps * std::pow(10, precision);
		int result = std::ceil(value * pow(10, precision));
		result = result - (rounder - result % rounder);
		value = result / (double)std::pow(10, precision);
	};
	auto roundup = [&] (double &value) {
		int rounder = info::mesh->contacts->eps * std::pow(10, precision);
		int result = std::floor(value * pow(10, precision));
		result = result + (rounder - result % rounder);
		value = result / (double)std::pow(10, precision);
	};

	rounddown(info::mesh->contacts->boundingBox[0].x);
	rounddown(info::mesh->contacts->boundingBox[0].y);
	rounddown(info::mesh->contacts->boundingBox[0].z);
	roundup(info::mesh->contacts->boundingBox[1].x);
	roundup(info::mesh->contacts->boundingBox[1].y);
	roundup(info::mesh->contacts->boundingBox[1].z);

	size_t xmaxsize = std::ceil((info::mesh->contacts->boundingBox[1].x - info::mesh->contacts->boundingBox[0].x) / info::mesh->contacts->eps);
	size_t ymaxsize = std::ceil((info::mesh->contacts->boundingBox[1].y - info::mesh->contacts->boundingBox[0].y) / info::mesh->contacts->eps);
	size_t zmaxsize = std::ceil((info::mesh->contacts->boundingBox[1].z - info::mesh->contacts->boundingBox[0].z) / info::mesh->contacts->eps);

	double avgsize = (xmaxsize + ymaxsize + zmaxsize) / 3.;
	info::mesh->contacts->groupsize = std::max(avgsize / std::pow(info::mesh->contacts->surface->enodes->structures(), 1. / 3) / 1.5, 1.);

	Communication::allReduce(&info::mesh->contacts->boundingBox[0], &info::mesh->contacts->globalBox[0], 3, MPI_DOUBLE, MPI_MIN);
	Communication::allReduce(&info::mesh->contacts->boundingBox[1], &info::mesh->contacts->globalBox[1], 3, MPI_DOUBLE, MPI_MAX);

	std::vector<std::vector<std::pair<esint, esint> > > grids(threads);

	double boxsize = info::mesh->contacts->eps * info::mesh->contacts->groupsize;
	size_t xsize = std::ceil((info::mesh->contacts->boundingBox[1].x - info::mesh->contacts->boundingBox[0].x) / boxsize);
	size_t ysize = std::ceil((info::mesh->contacts->boundingBox[1].y - info::mesh->contacts->boundingBox[0].y) / boxsize);

	// 2. map into grid
	///////////////////
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<esint, esint> > grid;

		int xmin, xmax, ymin, ymax, zmin, zmax;
		size_t gsize = 0;

		auto insert = [&] (Point &min, Point &max, esint eindex) {
			xmin = std::floor((min.x - info::mesh->contacts->boundingBox[0].x - epsilon) / boxsize);
			xmax = std::ceil((max.x - info::mesh->contacts->boundingBox[0].x + epsilon) / boxsize);
			ymin = std::floor((min.y - info::mesh->contacts->boundingBox[0].y - epsilon) / boxsize);
			ymax = std::ceil((max.y - info::mesh->contacts->boundingBox[0].y + epsilon) / boxsize);
			zmin = std::floor((min.z - info::mesh->contacts->boundingBox[0].z - epsilon) / boxsize);
			zmax = std::ceil((max.z - info::mesh->contacts->boundingBox[0].z + epsilon) / boxsize);

			grid.resize(gsize + (xmax - xmin) * (ymax - ymin) * (zmax - zmin), {0, eindex});
			for (int z = zmin; z < zmax; ++z) {
				for (int y = ymin; y < ymax; ++y) {
					for (int x = xmin; x < xmax; ++x, ++gsize) {
						grid[gsize].first = xsize * ysize * z + xsize * y + x;
					}
				}
			}
		};

		Point min, max;
		esint eindex = info::mesh->contacts->surface->edistribution[t];
		for (auto e = info::mesh->contacts->elements->cbegin(t); e != info::mesh->contacts->elements->cend(t); ++e, ++eindex) {
			min = max = e->front();
			for (auto n = e->begin() + 1; n != e->end(); ++n) {
				min.x = std::min(min.x, n->x);
				min.y = std::min(min.y, n->y);
				min.z = std::min(min.z, n->z);
				max.x = std::max(max.x, n->x);
				max.y = std::max(max.y, n->y);
				max.z = std::max(max.z, n->z);
			}
			insert(min, max, eindex);
		}

		grids[t].swap(grid);
	}

	{
		std::vector<size_t> gdistribution = { 0, grids[0].size() };
		for (size_t t = 1; t < threads; t++) {
			grids[0].insert(grids[0].end(), grids[t].begin(), grids[t].end());
			gdistribution.push_back(grids[0].size());
		}
		utils::sortWithInplaceMerge(grids[0], gdistribution);
	}

	// 3. create sparse structure
	///////////////////

	std::vector<size_t> distribution = tarray<size_t>::distribute(threads, grids[0].size());

	for (size_t t = 1; t < threads; t++) {
		while (distribution[t] != distribution[t + 1] && grids[0][distribution[t]].first == grids[0][distribution[t] - 1].first) {
			++distribution[t];
		}
	}

	std::vector<std::vector<esint> > gdist(threads), gdata(threads), gfilled(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> dist, data, filled;
		if (t == 0) {
			dist.push_back(0);
		}
		if (distribution[t] != distribution[t + 1]) {
			filled.push_back(grids[0][distribution[t]].first);
		}

		for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
			if (filled.back() != grids[0][i].first) {
				filled.push_back(grids[0][i].first);
				dist.push_back(data.size());
			}
			data.push_back(grids[0][i].second);
		}
		dist.push_back(data.size());

		gdist[t].swap(dist);
		gdata[t].swap(data);
		gfilled[t].swap(filled);
	}

	utils::threadDistributionToFullDistribution(gdist);
	for (size_t t = 1; t < threads; t++) {
		gdist[0].insert(gdist[0].end(), gdist[t].begin(), gdist[t].end());
		gdata[0].insert(gdata[0].end(), gdata[t].begin(), gdata[t].end());
		gfilled[0].insert(gfilled[0].end(), gfilled[t].begin(), gfilled[t].end());
	}

	// never threaded
	gdist.resize(1);
	gdata.resize(1);

	info::mesh->contacts->filledCells.swap(gfilled[0]);
	info::mesh->contacts->grid = new serializededata<esint, esint>(gdist, gdata);

	eslog::checkpointln("MESH: SURFACE LOCATIONS COMPUTED");
}

void computeContactNormals()
{
	if (info::mesh->contacts->elements == NULL) {
		computeSurfaceLocations();
	}

	if (info::mesh->contacts->enormals != NULL) {
		delete info::mesh->contacts->enormals;
	}
	info::mesh->contacts->enormals = new serializededata<esint, Point>(*info::mesh->contacts->elements);

//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const auto &epointers = info::mesh->surface->epointers->datatarray();
//		auto element = info::mesh->contacts->elements->begin(t);
//		auto normals = info::mesh->contacts->enormals->begin(t);
//		MatrixDense tangents(2, 3), coords(8, 3); // maximal size
//
//		for (auto e = info::mesh->contacts->surface->edistribution[t]; e < info::mesh->contacts->surface->edistribution[t + 1]; ++e, ++element, ++normals) {
//			coords.resize(epointers[e]->nodes, 3);
//			int cindex = 0;
//			for (auto n = element->begin(); n != element->end(); ++n, ++cindex) {
//				coords(cindex, 0) = n->x;
//				coords(cindex, 1) = n->y;
//				coords(cindex, 2) = n->z;
//			}
//			for (int n = 0; n < epointers[e]->nodes; n++) {
//				tangents.multiply((*epointers[e]->ndN)[n], coords);
//				Point t(tangents(0, 0), tangents(0, 1), tangents(0, 2));
//				Point s(tangents(1, 0), tangents(1, 1), tangents(1, 2));
//				normals->at(n) = Point::cross(t, s).normalize();
//			}
//		}
//	}

	std::vector<std::vector<Point> > sBuffer(info::mesh->contacts->gneighbors.size()), rBuffer(info::mesh->contacts->gneighbors.size());

	#pragma omp parallel for
	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		for (size_t i = 0; i < info::mesh->contacts->gnsurface[n].size(); ++i) {
			auto element = info::mesh->contacts->enormals->begin() + info::mesh->contacts->gnsurface[n][i];
			sBuffer[n].insert(sBuffer[n].end(), element->begin(), element->end());
		}
		rBuffer[n].resize(info::mesh->contacts->gnecoords[n]->datatarray().size());
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: exchange normals of geometric neighbors elements");
	}

	info::mesh->contacts->gnenormals.resize(info::mesh->contacts->gneighbors.size());
	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		info::mesh->contacts->gnenormals[n] = new serializededata<esint, Point>(tarray<esint>(info::mesh->contacts->gnecoords[n]->boundarytarray()), rBuffer[n]);
	}

	eslog::checkpointln("MESH: CONTACT NORMALS COMPUTED");
}

struct __contact__pair__ {
	// master, neigh, slave, start of triangles, number of triangles
	esint m, n, s, start, triangles;

	__contact__pair__(esint m, esint n, esint s, esint start): m(m), n(n), s(s), start(start), triangles(0) {}
};

void fireNormals()
{
	double threshold = 1e-6;
	Communication::serialize([&] () {
		return;
		std::cout << "\n";
		std::cout << " -- " << info::mpi::rank << " -- \n";
//		std::cout << "etypes: ";
//		for (auto etype = info::mesh->contacts->epointers->datatarray().begin(); etype != info::mesh->contacts->epointers->datatarray().end(); ++etype) {
//			std::cout << (int)(*etype)->code << " ";
//		}
		std::cout << "elements: " << *info::mesh->contacts->elements << "\n";
		std::cout << "normals: " << *info::mesh->contacts->enormals << "\n";
		std::cout << "my close: " << *info::mesh->contacts->closeElements << "\n";
		//		std::cout << "neighs: " << *info::mesh->contacts->surface->neighbors << "\n";
		for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
			std::cout << " --- with " << info::mesh->contacts->gneighbors[n] << " ---\n";
//			std::cout << "etypes [" << info::mesh->contacts->gneighbors[n] << "]: ";
//			for (auto etype = info::mesh->contacts->gnepointers[n]->datatarray().begin(); etype != info::mesh->contacts->gnepointers[n]->datatarray().end(); ++etype) {
//				std::cout << (int)(*etype)->code << " ";
//			}
			std::cout << "elements: " <<* info::mesh->contacts->gnecoords[n] << "\n";
			std::cout << "normals: " <<* info::mesh->contacts->gnenormals[n] << "\n";
			std::cout << "IDs: " << *info::mesh->contacts->gneIDs[n] << "\n";
//			std::cout << "grid: " << *info::mesh->contacts->gngrid[n] << "\n";
			std::cout << "surface: " << info::mesh->contacts->gnsurface[n];
			std::cout << "close: " << *info::mesh->contacts->gncloseElements[n] << "\n";
		}
	});

	std::vector<Point> p1, p2;
	std::vector<std::vector<Point> > intersection;

	esint clips = 0;

	double time = 0;
	int counter = 0;

	std::vector<__contact__pair__> pairs;
	std::vector<Triangle> triangles;

	auto neighbors = info::mesh->contacts->surface->neighbors->begin(); // not used now
	auto epoints = info::mesh->contacts->elements->begin();
	auto closest = info::mesh->contacts->closeElements->begin();
	auto enormals = info::mesh->contacts->enormals->begin();
	for (size_t e = 0; e < info::mesh->contacts->surface->enodes->structures(); ++e, ++epoints, ++neighbors, ++closest, ++enormals) {
		Point &normal = enormals->at(0);
		Point center;
		for (auto p = epoints->begin(); p != epoints->end(); ++p) {
			center += *p;
		}
		center /= epoints->size();
		bool max_x = false, max_y = false, max_z = false;
		if (std::fabs(normal.x) <= std::fabs(normal.z) && std::fabs(normal.y) <= std::fabs(normal.z)) {
			max_z = true; max_y = false; max_x = false;
		}
		if (std::fabs(normal.x) <= std::fabs(normal.y) && std::fabs(normal.z) <= std::fabs(normal.y)) {
			max_y = true; max_x = false; max_z = false;
		}
		if (std::fabs(normal.y) <= std::fabs(normal.x) && std::fabs(normal.z) <= std::fabs(normal.x)) {
			max_x = true; max_y = false; max_z = false;
		}
		p1.clear();
		for (auto p = epoints->begin(); p != epoints->end(); ++p) {
			Point pp = *p - normal * (*p * normal);
			if (max_x) { p1.push_back(Point(pp.y, pp.z, 0)); }
			if (max_y) { p1.push_back(Point(pp.x, pp.z, 0)); }
			if (max_z) { p1.push_back(Point(pp.x, pp.y, 0)); }
		}
		for (auto c = closest->begin(); c != closest->end(); ++c) {
			if ((esint)e < *c) {
				auto other = info::mesh->contacts->elements->begin() + *c;
				Point &onormal = (info::mesh->contacts->enormals->begin() + *c)->at(0);
				Point ocenter;
				for (auto p = other->begin(); p != other->end(); ++p) {
					ocenter += *p;
				}
				ocenter /= other->size();
				double distance = (ocenter - center) * normal;
				if (distance < -0.1 || distance > 0.5 || onormal * normal > threshold || std::fabs(onormal * normal) < threshold) {
					continue;
				} else {
					p2.clear();
					for (auto p = other->begin(); p != other->end(); ++p) {
						Point pp = *p - normal * (*p * normal);
						if (max_x) { p2.push_back(Point(pp.y, pp.z, 0)); }
						if (max_y) { p2.push_back(Point(pp.x, pp.z, 0)); }
						if (max_z) { p2.push_back(Point(pp.x, pp.y, 0)); }
					}
					Point pcenter;
					for (size_t pp = 0; pp < p1.size(); ++pp) {
						pcenter += p1[pp];
					}
					pcenter /= p1.size();
					if (Point::cross(Point(p1[0] - pcenter), Point(p1[1] - pcenter)).z < 0) {
						p1 = std::vector<Point>(p1.rbegin(), p1.rend());
					}
					pcenter = Point();
					for (size_t pp = 0; pp < p2.size(); ++pp) {
						pcenter += p2[pp];
					}
					pcenter /= p2.size();
					if (Point::cross(Point(p2[0] - pcenter), Point(p2[1] - pcenter)).z < 0) {
						p2 = std::vector<Point>(p2.rbegin(), p2.rend());
					}

					intersection.clear();
					double start = TimeLogger::time();
//					dummyClip(p1, p2, intersection);
//					if (info::mpi::rank) {
//						clip(p1, p2, intersection);
//					} else {
//						dummyClip(p1, p2, intersection);
//					}
					time += TimeLogger::time() - start;
					++counter;
					if (intersection.size()) {
						++clips;
						pairs.push_back(__contact__pair__( e, -1, *c, triangles.size()));
						for (size_t i = 0; i < intersection.size(); i++) {
							pairs.back().triangles += intersection[i].size() - 2;
							for (size_t j = 0; j < intersection[i].size(); j++) {
								if (max_x) { intersection[i][j] = Point(epoints->front().x, intersection[i][j].x, intersection[i][j].y); }
								if (max_y) { intersection[i][j] = Point(intersection[i][j].x, epoints->front().y, intersection[i][j].y); }
								if (max_z) { intersection[i][j] = Point(intersection[i][j].x, intersection[i][j].y, epoints->front().z); }
							}
							for (size_t j = 2; j < intersection[i].size(); j++) {
								triangles.push_back(Triangle{ intersection[i][0], intersection[i][j - 1], intersection[i][j] });
							}
						}
					}
				}
			}
		}
	}

	{ // NEIGHBORING CONTACTS
		info::mesh->contacts->gnintersection.resize(info::mesh->contacts->gneighbors.size());
		for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) { // neighboring elements
			auto nclosest = info::mesh->contacts->gncloseElements[n]->cbegin();
			for (size_t i = 0; i < info::mesh->contacts->gnsurface[n].size(); ++i, ++nclosest) {
				auto mypoints = info::mesh->contacts->elements->begin() + info::mesh->contacts->gnsurface[n][i];
				auto mynormals = info::mesh->contacts->enormals->begin() + info::mesh->contacts->gnsurface[n][i];
				Point &normal = mynormals->at(0);
				Point center;
				for (auto p = mypoints->begin(); p != mypoints->end(); ++p) {
					center += *p;
				}
				center /= mypoints->size();

				bool max_x = false, max_y = false, max_z = false;
				if (std::fabs(normal.x) <= std::fabs(normal.z) && std::fabs(normal.y) <= std::fabs(normal.z)) {
					max_z = true; max_y = false; max_x = false;
				}
				if (std::fabs(normal.x) <= std::fabs(normal.y) && std::fabs(normal.z) <= std::fabs(normal.y)) {
					max_y = true; max_x = false; max_z = false;
				}
				if (std::fabs(normal.y) <= std::fabs(normal.x) && std::fabs(normal.z) <= std::fabs(normal.x)) {
					max_x = true; max_y = false; max_z = false;
				}

				p1.clear();
				for (auto p = mypoints->begin(); p != mypoints->end(); ++p) {
					Point pp = *p - normal * (*p * normal);
					if (max_x) { p1.push_back(Point(pp.y, pp.z, 0)); }
					if (max_y) { p1.push_back(Point(pp.x, pp.z, 0)); }
					if (max_z) { p1.push_back(Point(pp.x, pp.y, 0)); }
				}

				for (auto ne = nclosest->begin(); ne != nclosest->end(); ++ne) {
					auto other = info::mesh->contacts->gnecoords[n]->begin() + *ne;
					Point &onormal = (info::mesh->contacts->gnenormals[n]->begin() + *ne)->at(0);
					Point ocenter;
					for (auto p = other->begin(); p != other->end(); ++p) {
						ocenter += *p;
					}
					ocenter /= other->size();
					double distance = (ocenter - center) * normal;
					if (distance < -0.1 || distance > 0.5 || onormal * normal > threshold || std::fabs(onormal * normal) < threshold) {
						continue;
					} else {
						p2.clear();
						for (auto p = other->begin(); p != other->end(); ++p) {
							Point pp = *p - normal * (*p * normal);
							if (max_x) { p2.push_back(Point(pp.y, pp.z, 0)); }
							if (max_y) { p2.push_back(Point(pp.x, pp.z, 0)); }
							if (max_z) { p2.push_back(Point(pp.x, pp.y, 0)); }
						}
						Point pcenter;
						for (size_t pp = 0; pp < p1.size(); ++pp) {
							pcenter += p1[pp];
						}
						pcenter /= p1.size();
						if (Point::cross(Point(p1[0] - pcenter), Point(p1[1] - pcenter)).z < 0) {
							p1 = std::vector<Point>(p1.rbegin(), p1.rend());
						}
						pcenter = Point();
						for (size_t pp = 0; pp < p2.size(); ++pp) {
							pcenter += p2[pp];
						}
						pcenter /= p2.size();
						if (Point::cross(Point(p2[0] - pcenter), Point(p2[1] - pcenter)).z < 0) {
							p2 = std::vector<Point>(p2.rbegin(), p2.rend());
						}

						intersection.clear();
						double start = TimeLogger::time();
//						dummyClip(p1, p2, intersection);
//						if (info::mpi::rank) {
//							clip(p1, p2, intersection);
//						} else {
//							dummyClip(p1, p2, intersection);
//						}
						time += TimeLogger::time() - start;
						++counter;

						if (intersection.size()) {
							++clips;
							pairs.push_back(__contact__pair__(info::mesh->contacts->gnsurface[n][i], n, *ne, triangles.size()));
							for (size_t i = 0; i < intersection.size(); i++) {
								pairs.back().triangles += intersection[i].size() - 2;
								for (size_t j = 0; j < intersection[i].size(); j++) {
									if (max_x) { intersection[i][j] = Point(mypoints->front().x, intersection[i][j].x, intersection[i][j].y); }
									if (max_y) { intersection[i][j] = Point(intersection[i][j].x, mypoints->front().y, intersection[i][j].y); }
									if (max_z) { intersection[i][j] = Point(intersection[i][j].x, intersection[i][j].y, mypoints->front().z); }
								}
								for (size_t j = 2; j < intersection[i].size(); j++) {
									triangles.push_back(Triangle{ intersection[i][0], intersection[i][j - 1], intersection[i][j] });
								}
							}
						}
					}
				}
			}
		}
	}

	if (pairs.size()) {
		std::vector<esint> permutation(pairs.size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
			return pairs[i].m < pairs[j].m;
		});

		std::vector<esint> pairsdist, pairsdata, tdist(1);
		std::vector<Triangle> tdata;

		esint prev = -1;
		for (auto i = permutation.begin(); i != permutation.end(); ++i) {
			if (pairs[*i].m != prev) {
				pairsdist.push_back(pairsdata.size());
				pairsdata.push_back(prev = pairs[*i].m);
				pairsdata.push_back(0);
			}
			++pairsdata[pairsdist.back() + 1];
			pairsdata.push_back(pairs[*i].n);
			pairsdata.push_back(pairs[*i].s);
			for (esint t = 0; t < pairs[*i].triangles; ++t) {
				tdata.push_back(triangles[pairs[*i].start + t]);
			}
			tdist.push_back(tdata.size());
		}
		pairsdist.push_back(pairsdata.size());

		info::mesh->contacts->contactPairs = new serializededata<esint, esint>(pairsdist, pairsdata);
		info::mesh->contacts->intersections = new serializededata<esint, Triangle>(tdist, tdata);
	}

//	Communication::serialize([&] () {
//		std::cout << " -- " << info::mpi::rank << " -- \n";
//		for (size_t i = 0; i < pairs.size(); ++i) {
//			std::cout << pairs[i].m << " -> " << pairs[i].s << "[" << pairs[i].triangles << "]\n";
//		}
//		std::cout << *info::mesh->contacts->contactPairs << "\n";
//		std::cout << *info::mesh->contacts->intersections << "\n";
//	});


	eslog::checkpointln("MESH: FIRE NORMALS");

	DebugOutput::contact();
}


void findCloseElements()
{
	size_t threads = info::env::OMP_NUM_THREADS;
	double epsilon = 1e-6;

	double boxsize = info::mesh->contacts->eps * info::mesh->contacts->groupsize;
	int xsize = std::ceil((info::mesh->contacts->boundingBox[1].x - info::mesh->contacts->boundingBox[0].x) / boxsize);
	int ysize = std::ceil((info::mesh->contacts->boundingBox[1].y - info::mesh->contacts->boundingBox[0].y) / boxsize);
	int zsize = std::ceil((info::mesh->contacts->boundingBox[1].z - info::mesh->contacts->boundingBox[0].z) / boxsize);

	auto areIntersected = [&] (Point *block1, Point *block2) {
		return	!
				(block1[1].x + epsilon < block2[0].x || block2[1].x + epsilon < block1[0].x) ||
				(block1[1].y + epsilon < block2[0].y || block2[1].y + epsilon < block1[0].y) ||
				(block1[1].z + epsilon < block2[0].z || block2[1].z + epsilon < block1[0].z);
	};

	auto intersect = [&] (Point *intersection, Point *block1, Point *block2) {
		intersection[0].x = std::max(block1[0].x, block2[0].x);
		intersection[0].y = std::max(block1[0].y, block2[0].y);
		intersection[0].z = std::max(block1[0].z, block2[0].z);

		intersection[1].x = std::min(block1[1].x, block2[1].x);
		intersection[1].y = std::min(block1[1].y, block2[1].y);
		intersection[1].z = std::min(block1[1].z, block2[1].z);
	};

	std::vector<Point> boxes(2 * info::mpi::size);

	Communication::allGather(info::mesh->contacts->boundingBox, boxes.data(), 6, MPI_DOUBLE);

	for (int r = 0; r < info::mpi::size; r++) {
		if (r != info::mpi::rank) {
			if (areIntersected(info::mesh->contacts->boundingBox, boxes.data() + 2 * r)) {
				info::mesh->contacts->gneighbors.push_back(r);
			}
		}
	}

	// EXCHANGE GROUP SIZES
	std::vector<std::vector<size_t> > sGroupSize(info::mesh->contacts->gneighbors.size(), { info::mesh->contacts->groupsize }), rGroupSize(info::mesh->contacts->gneighbors.size(), { 0 });
	if (!Communication::exchangeKnownSize(sGroupSize, rGroupSize, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: exchange group size.\n");
	}

	std::vector<std::vector<esint> > sFilled(info::mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sBlock(info::mesh->contacts->gneighbors.size()), rBlock(info::mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sOffsets(info::mesh->contacts->gneighbors.size()), rOffsets(info::mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sIDs(info::mesh->contacts->gneighbors.size()), rIDs(info::mesh->contacts->gneighbors.size());
	std::vector<std::vector<int> > sTypes(info::mesh->contacts->gneighbors.size()), rTypes(info::mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sDist(info::mesh->contacts->gneighbors.size()), rDist(info::mesh->contacts->gneighbors.size());
	std::vector<std::vector<Point> > sData(info::mesh->contacts->gneighbors.size()), rData(info::mesh->contacts->gneighbors.size());
	std::vector<std::vector<esint> > sNodeData(info::mesh->contacts->gneighbors.size()), rNodeData(info::mesh->contacts->gneighbors.size());

	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		Point intersection[2];
		intersect(intersection, info::mesh->contacts->boundingBox, boxes.data() + 2 * info::mesh->contacts->gneighbors[n]);

		int xmin = std::max((int)std::floor((intersection[0].x - info::mesh->contacts->boundingBox[0].x - epsilon) / boxsize), 0);
		int ymin = std::max((int)std::floor((intersection[0].y - info::mesh->contacts->boundingBox[0].y - epsilon) / boxsize), 0);
		int zmin = std::max((int)std::floor((intersection[0].z - info::mesh->contacts->boundingBox[0].z - epsilon) / boxsize), 0);
		int xmax = std::min((int)std::ceil ((intersection[1].x - info::mesh->contacts->boundingBox[0].x + epsilon) / boxsize), xsize);
		int ymax = std::min((int)std::ceil ((intersection[1].y - info::mesh->contacts->boundingBox[0].y + epsilon) / boxsize), xsize);
		int zmax = std::min((int)std::ceil ((intersection[1].z - info::mesh->contacts->boundingBox[0].z + epsilon) / boxsize), xsize);

		auto begin = std::lower_bound(info::mesh->contacts->filledCells.begin(), info::mesh->contacts->filledCells.end(), xsize * ysize * zmin);
		auto end   = std::lower_bound(info::mesh->contacts->filledCells.begin(), info::mesh->contacts->filledCells.end(), xsize * ysize * zmax);
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, end - begin);
		for (size_t t = 0; t <= threads; t++) {
			distribution[t] += begin - info::mesh->contacts->filledCells.begin();
		}

		std::vector<std::vector<esint> > filled(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> tfilled;
			esint y, x;
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
				y = info::mesh->contacts->filledCells[i] % (xsize * ysize) / xsize;
				x = info::mesh->contacts->filledCells[i] % xsize;
				if (xmin <= x && x < xmax && ymin <= y && y < ymax) {
					tfilled.push_back(i);
				}
			}

			filled[t].swap(tfilled);
		}

		sFilled[n].clear();
		for (size_t t = 0; t < threads; t++) {
			sFilled[n].insert(sFilled[n].end(), filled[t].begin(), filled[t].end());
		}
	}

	info::mesh->contacts->gnsurface.resize(info::mesh->contacts->gneighbors.size());

	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, sFilled[n].size());
		sBlock[n].resize(sFilled[n].size() + 1);

		std::vector<std::vector<esint> > tIDs(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> IDs;

			if (distribution[t] < distribution[t + 1]) {
				esint prevBlock = sFilled[n][distribution[t]];
				auto cell = info::mesh->contacts->grid->cbegin() + prevBlock;
				for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
					cell += sFilled[n][i] - prevBlock;
					sFilled[n][i] = info::mesh->contacts->filledCells[prevBlock = sFilled[n][i]];
					IDs.insert(IDs.end(), cell->begin(), cell->end());
					sBlock[n][i + 1] = IDs.size();
				}

				tIDs[t].swap(IDs);
			}
		}

		utils::threadDistributionToFullDistribution(sBlock[n], distribution);
		for (size_t t = 0; t < threads; t++) {
			distribution[t] = sOffsets[n].size();
			sOffsets[n].insert(sOffsets[n].end(), tIDs[t].begin(), tIDs[t].end());
		}
		distribution[threads] = sOffsets[n].size();
		sIDs[n] = sOffsets[n];

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			utils::sortAndRemoveDuplicates(tIDs[t]);
		}
		utils::sortAndRemoveDuplicates(sIDs[n]);

		for (size_t t = 1; t < threads; t++) {
			tIDs[0].insert(tIDs[0].end(), tIDs[t].begin(), tIDs[t].end());
		}
		utils::sortAndRemoveDuplicates(tIDs[0]);
		info::mesh->contacts->gnsurface[n].swap(tIDs[0]);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
				sOffsets[n][i] = std::lower_bound(info::mesh->contacts->gnsurface[n].begin(), info::mesh->contacts->gnsurface[n].end(), sOffsets[n][i]) - info::mesh->contacts->gnsurface[n].begin();
			}
		}
		for (size_t i = 0; i < sIDs[n].size(); ++i) {
			sIDs[n][i] = info::mesh->contacts->surface->IDs->datatarray()[sIDs[n][i]];
		}

		distribution = tarray<size_t>::distribute(threads, info::mesh->contacts->gnsurface[n].size());

		std::vector<std::vector<esint> > dist(threads), types(threads), nodes(threads);
		std::vector<std::vector<Point> > data(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> tdist, ttypes, tnodes;
			std::vector<Point> tdata;
			if (t == 0) {
				tdist.reserve(distribution[t + 1] - distribution[t] + 1);
				tdist.push_back(0);
			} else {
				tdist.reserve(distribution[t + 1] - distribution[t]);
			}
			if (distribution[t] < distribution[t + 1]) {
				esint prev = info::mesh->contacts->gnsurface[n][distribution[t]];
				auto face = info::mesh->contacts->elements->cbegin() + prev;
				auto fnodes = info::mesh->contacts->enodes->cbegin() + prev;
				for (size_t i = distribution[t]; i < distribution[t + 1]; prev = info::mesh->contacts->gnsurface[n][i++]) {
					face += info::mesh->contacts->gnsurface[n][i] - prev;
					fnodes += info::mesh->contacts->gnsurface[n][i] - prev;
					tdata.insert(tdata.end(), face->begin(), face->end());
					tnodes.insert(tnodes.end(), fnodes->begin(), fnodes->end());
					tdist.push_back(tdata.size());
					ttypes.push_back((int)info::mesh->contacts->epointers->datatarray()[info::mesh->contacts->gnsurface[n][i]]->code);
				}
			}
			dist[t].swap(tdist);
			data[t].swap(tdata);
			nodes[t].swap(tnodes);
			types[t].swap(ttypes);
		}
		utils::threadDistributionToFullDistribution(dist);

		for (size_t t = 0; t < threads; t++) {
			sDist[n].insert(sDist[n].end(), dist[t].begin(), dist[t].end());
			sData[n].insert(sData[n].end(), data[t].begin(), data[t].end());
			sNodeData[n].insert(sNodeData[n].end(), nodes[t].begin(), nodes[t].end());
			sTypes[n].insert(sTypes[n].end(), types[t].begin(), types[t].end());
		}
	}

	info::mesh->contacts->gnfilled.resize(info::mesh->contacts->gneighbors.size());
	info::mesh->contacts->gneIDs.resize(info::mesh->contacts->gneighbors.size());
	info::mesh->contacts->gngrid.resize(info::mesh->contacts->gneighbors.size());
	info::mesh->contacts->gnecoords.resize(info::mesh->contacts->gneighbors.size());
	info::mesh->contacts->gnenodes.resize(info::mesh->contacts->gneighbors.size());
	info::mesh->contacts->gnepointers.resize(info::mesh->contacts->gneighbors.size());

	if (!Communication::exchangeUnknownSize(sFilled, info::mesh->contacts->gnfilled, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - filled boxed.\n");
	}
	if (!Communication::exchangeUnknownSize(sBlock, rBlock, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - blocks sizes.\n");
	}
	if (!Communication::exchangeUnknownSize(sIDs, rIDs, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - IDs.\n");
	}
	if (!Communication::exchangeUnknownSize(sOffsets, rOffsets, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - Offsets.\n");
	}
	if (!Communication::exchangeUnknownSize(sDist, rDist, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - faces distribution.\n");
	}
	if (!Communication::exchangeUnknownSize(sData, rData, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - faces data.\n");
	}
	if (!Communication::exchangeUnknownSize(sNodeData, rNodeData, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - faces node data.\n");
	}
	if (!Communication::exchangeUnknownSize(sTypes, rTypes, info::mesh->contacts->gneighbors)) {
		eslog::error("ESPRESO internal error: contacts - etypoes.\n");
	}

	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); ++n) {
		info::mesh->contacts->gneIDs[n] = new serializededata<esint, esint>(1, rIDs[n]);
		info::mesh->contacts->gngrid[n] = new serializededata<esint, esint>(rBlock[n], rOffsets[n]);
		info::mesh->contacts->gnecoords[n] = new serializededata<esint, Point>(rDist[n], rData[n]);
		info::mesh->contacts->gnenodes[n] = new serializededata<esint, esint>(rDist[n], rNodeData[n]);

		std::vector<Element*> epointers;
		for (size_t i = 0; i < rTypes[n].size(); i++) {
			epointers.push_back(&Mesh::edata[rTypes[n][i]]);
		}
		info::mesh->contacts->gnepointers[n] = new serializededata<esint, Element*>(1, epointers);
	}

	std::vector<std::vector<esint> > closeElementsDist(threads), closeElementsData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tcloseElementsDist, tcloseElementsData;
		if (t == 0) {
			tcloseElementsDist.reserve(info::mesh->contacts->elements->cbegin(t)->begin() - info::mesh->contacts->elements->cbegin(t)->begin() + 1);
			tcloseElementsDist.push_back(0);
		} else {
			tcloseElementsDist.reserve(info::mesh->contacts->elements->cbegin(t)->begin() - info::mesh->contacts->elements->cbegin(t)->begin());
		}
		Point min, max;
		int xmin, xmax, ymin, ymax, zmin, zmax;
		size_t prevsize;
		std::vector<esint>::const_iterator zbegin, zend, ybegin, yend, xbegin, xend;
		for (auto e = info::mesh->contacts->elements->cbegin(t); e != info::mesh->contacts->elements->cend(t); ++e) {
			min = max = e->front();
			for (auto n = e->begin() + 1; n != e->end(); ++n) {
				min.x = std::min(min.x, n->x);
				min.y = std::min(min.y, n->y);
				min.z = std::min(min.z, n->z);
				max.x = std::max(max.x, n->x);
				max.y = std::max(max.y, n->y);
				max.z = std::max(max.z, n->z);
			}

			xmin = std::max((int)std::floor((min.x - info::mesh->contacts->boundingBox[0].x - info::mesh->contacts->eps - epsilon) / boxsize), 0);
			ymin = std::max((int)std::floor((min.y - info::mesh->contacts->boundingBox[0].y - info::mesh->contacts->eps - epsilon) / boxsize), 0);
			zmin = std::max((int)std::floor((min.z - info::mesh->contacts->boundingBox[0].z - info::mesh->contacts->eps - epsilon) / boxsize), 0);
			xmax = std::min((int)std::ceil ((max.x - info::mesh->contacts->boundingBox[0].x + info::mesh->contacts->eps + epsilon) / boxsize), xsize);
			ymax = std::min((int)std::ceil ((max.y - info::mesh->contacts->boundingBox[0].y + info::mesh->contacts->eps + epsilon) / boxsize), ysize);
			zmax = std::min((int)std::ceil ((max.z - info::mesh->contacts->boundingBox[0].z + info::mesh->contacts->eps + epsilon) / boxsize), zsize);

			auto begin = info::mesh->contacts->filledCells.begin();
			auto end   = std::lower_bound(info::mesh->contacts->filledCells.begin(), info::mesh->contacts->filledCells.end(), xsize * ysize * (zmax - 1) + xsize * (ymax - 1) + (xmax - 1) + 1);
			auto faces = info::mesh->contacts->grid->cbegin();

			auto prevcell = begin;
			prevsize = tcloseElementsData.size();
			for (int z = zmin; begin != end && z < zmax; ++z) {
				for (int y = ymin; begin != end && y < ymax; ++y) {
					begin = std::lower_bound(begin, end, xsize * ysize * z + xsize * y + xmin);
					if (begin != end) {
						faces += begin - prevcell;
						prevcell = begin;
					} else {
						break;
					}
					while (*begin < xsize * ysize * z + xsize * y + xmax) {
						tcloseElementsData.insert(tcloseElementsData.end(), faces->begin(), faces->end());
						++begin;
						if (begin != end) {
							faces += begin - prevcell;
							prevcell = begin;
						} else {
							break;
						}
					}
				}
			}
			utils::sortAndRemoveDuplicates(tcloseElementsData, prevsize);
			tcloseElementsDist.push_back(tcloseElementsData.size());
		}

		closeElementsDist[t].swap(tcloseElementsDist);
		closeElementsData[t].swap(tcloseElementsData);
	}

	utils::threadDistributionToFullDistribution(closeElementsDist);

	info::mesh->contacts->closeElements = new serializededata<esint, esint>(closeElementsDist, closeElementsData);
	info::mesh->contacts->gncloseElements.resize(info::mesh->contacts->gneighbors.size());
	for (size_t n = 0; n < info::mesh->contacts->gneighbors.size(); n++) {
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, info::mesh->contacts->gnsurface[n].size());
		std::vector<std::vector<esint> > ncloseElementsDist(threads), ncloseElementsData(threads);

		double nboxsize = info::mesh->contacts->eps * rGroupSize[n].front();
		int nxsize = std::ceil((boxes[2 * info::mesh->contacts->gneighbors[n] + 1].x - boxes[2 * info::mesh->contacts->gneighbors[n]].x) / nboxsize);
		int nysize = std::ceil((boxes[2 * info::mesh->contacts->gneighbors[n] + 1].y - boxes[2 * info::mesh->contacts->gneighbors[n]].y) / nboxsize);
		int nzsize = std::ceil((boxes[2 * info::mesh->contacts->gneighbors[n] + 1].z - boxes[2 * info::mesh->contacts->gneighbors[n]].z) / nboxsize);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point min, max;
			int xmin, xmax, ymin, ymax, zmin, zmax;
			size_t prevsize;
			std::vector<esint> tcloseElementsDist, tcloseElementsData;
			if (t == 0) {
				tcloseElementsDist.reserve(distribution[t + 1] - distribution[t] + 1);
				tcloseElementsDist.push_back(0);
			} else {
				tcloseElementsDist.reserve(distribution[t + 1] - distribution[t]);
			}

			if (distribution[t] != distribution[t + 1]) {
				esint prev = info::mesh->contacts->gnsurface[n][distribution[t]];
				auto e = info::mesh->contacts->elements->cbegin() + prev;
				for (size_t i = distribution[t]; i < distribution[t + 1]; prev = info::mesh->contacts->gnsurface[n][i++]) {
					e += info::mesh->contacts->gnsurface[n][i] - prev;
					min = max = e->front();
					for (auto nn = e->begin() + 1; nn != e->end(); ++nn) {
						min.x = std::min(min.x, nn->x);
						min.y = std::min(min.y, nn->y);
						min.z = std::min(min.z, nn->z);
						max.x = std::max(max.x, nn->x);
						max.y = std::max(max.y, nn->y);
						max.z = std::max(max.z, nn->z);
					}

					xmin = std::max((int)std::floor((min.x - boxes[2 * info::mesh->contacts->gneighbors[n]].x - info::mesh->contacts->eps - epsilon) / nboxsize), 0);
					ymin = std::max((int)std::floor((min.y - boxes[2 * info::mesh->contacts->gneighbors[n]].y - info::mesh->contacts->eps - epsilon) / nboxsize), 0);
					zmin = std::max((int)std::floor((min.z - boxes[2 * info::mesh->contacts->gneighbors[n]].z - info::mesh->contacts->eps - epsilon) / nboxsize), 0);
					xmax = std::min((int)std::ceil ((max.x - boxes[2 * info::mesh->contacts->gneighbors[n]].x + info::mesh->contacts->eps + epsilon) / nboxsize), nxsize);
					ymax = std::min((int)std::ceil ((max.y - boxes[2 * info::mesh->contacts->gneighbors[n]].y + info::mesh->contacts->eps + epsilon) / nboxsize), nysize);
					zmax = std::min((int)std::ceil ((max.z - boxes[2 * info::mesh->contacts->gneighbors[n]].z + info::mesh->contacts->eps + epsilon) / nboxsize), nzsize);

					auto begin = info::mesh->contacts->gnfilled[n].begin();
					auto end   = std::lower_bound(info::mesh->contacts->gnfilled[n].begin(), info::mesh->contacts->gnfilled[n].end(), nxsize * nysize * (zmax - 1) + nxsize * (ymax - 1) + (xmax - 1) + 1);
					auto faces = info::mesh->contacts->gngrid[n]->cbegin();

					auto prevcell = begin;
					prevsize = tcloseElementsData.size();
					for (int z = zmin; begin != end && z < zmax; ++z) {
						for (int y = ymin; begin != end && y < ymax; ++y) {
							begin = std::lower_bound(begin, end, nxsize * nysize * z + nxsize * y + xmin);
							if (begin != end) {
								faces += begin - prevcell;
								prevcell = begin;
							} else {
								break;
							}
							while (*begin < nxsize * nysize * z + nxsize * y + xmax) {
								tcloseElementsData.insert(tcloseElementsData.end(), faces->begin(), faces->end());
								++begin;
								if (begin != end) {
									faces += begin - prevcell;
									prevcell = begin;
								} else {
									break;
								}
							}
						}
					}
					utils::sortAndRemoveDuplicates(tcloseElementsData, prevsize);
					tcloseElementsDist.push_back(tcloseElementsData.size());
				}
			}

			ncloseElementsDist[t].swap(tcloseElementsDist);
			ncloseElementsData[t].swap(tcloseElementsData);
		}

		utils::threadDistributionToFullDistribution(ncloseElementsDist);

		info::mesh->contacts->gncloseElements[n] = new serializededata<esint, esint>(ncloseElementsDist, ncloseElementsData);
	}

	DebugOutput::closeElements();

	eslog::checkpointln("MESH: CONTACT INTERFACE SEARCHED");
}

int intersection(Point &p0, Point &p1, Point &q0, Point &q1, Point &first, Point &second)
{
	double epsilon = 1e-8;
	auto perp = [] (Point &u, Point &v) {
		return u.x * v.y - u.y * v.x;
	};

	Point u = p1 - p0;
	Point v = q1 - q0;
	Point w = p0 - q0;
	double D = perp(u, v);

	// test if they are parallel (includes either being a point)
	if (fabs(D) < epsilon) { // lines are parallel
		if (perp(u, w) != 0 || perp(v, w) != 0) {
			return 0; // they are NOT collinear
		}
		// they are collinear

		// they are collinear segments - get  overlap (or not)
		double t0, t1; // endpoints of S1 in eqn for S2
		Point w2 = p1 - q0;
		if (v.x != 0) {
			t0 = w.x / v.x;
			t1 = w2.x / v.x;
		} else {
			t0 = w.y / v.y;
			t1 = w2.y / v.y;
		}
		if (t0 > t1) { // must have t0 smaller than t1
			double t = t0;
			t0 = t1;
			t1 = t;  // swap if not
		}
		if (t0 > 1 || t1 < 0) {
			return 0; // NO overlap
		}
		t0 = t0 < 0 ? 0 : t0; // clip to min 0
		t1 = t1 > 1 ? 1 : t1; // clip to max 1
		if (t0 == t1) { // intersect is a point
			first = q0 + v * t0;
			return 1;
		}

		// they overlap in a valid subsegment
		first = q0 + v * t0;
		first.z = t0;
		second = q0 + v * t1;
		second.z = t1;
		return 2;
	}

	// the segments are skew and may intersect in a point
	// get the intersect parameter for S1
	double s = perp(v, w) / D;
	if (s < 0 || s > 1) {
		return 0; // no intersect with S1
	}

	// get the intersect parameter for S2
	double t = perp(u, w) / D;
	if (t < 0 || t > 1) {
		return 0; // no intersect with S2
	}

	first = q0 + v * t; // compute S1 intersect point
	first.z = t;
	return 1;
}

int isIn(Point &p, std::vector<Point> &polygon)
{
	int n = polygon.size(), cn = 0;    // the crossing number counter

	for (size_t i = 0; i < polygon.size(); i++) {
		if (
				(polygon[i].y <= p.y && polygon[(i + 1) % n].y >  p.y) || // an upward crossing
				(polygon[i].y >  p.y && polygon[(i + 1) % n].y <= p.y)    // a downward crossing
			) {
			// compute the actual edge-ray intersect x-coordinate
			double vt = (p.y - polygon[i].y) / (polygon[(i + 1) % n].y - polygon[i].y);
			if (p.x < polygon[i].x + vt * (polygon[(i + 1) % n].x - polygon[i].x)) {
				++cn;   // a valid crossing of y = P.y right of P.x
			}
		}
	}
	return (cn & 1);    // 0 if even (out), and 1 if odd (in)
}

struct __line__ {
	Point start, end;
};

struct __line_iterator__ {
	std::vector<Point> &polygon;
	std::vector<Point>::iterator start;
	std::vector<Point>::iterator end;

	__line_iterator__(std::vector<Point> &polygon)
	: polygon(polygon), start(polygon.begin()), end(polygon.begin() + 1) {}

	void next() {
		if (++start == polygon.end()) { start = polygon.begin(); }
		if (++end == polygon.end()) { end = polygon.begin(); }
	}
};

void dummyClip(std::vector<Point> &p, std::vector<Point> &q, std::vector<std::vector<Point> > &res)
{
	auto addInnerLines = [&] (__line_iterator__ current, __line_iterator__ other, std::vector<__line__> &lines) {
		std::vector<Point> intersections;
		for (size_t i = 0; i < current.polygon.size(); ++i, current.next()) {
			intersections.clear();
			if (isIn(*current.start, other.polygon)) {
				intersections.push_back(*current.start);
				intersections.back().z = 0;
			}
			if (isIn(*current.end, other.polygon)) {
				intersections.push_back(*current.end);
				intersections.back().z = 1;
			}

			Point first, second;
			double last = -1;
			for (size_t j = 0; j < other.polygon.size(); j++, other.next()) {
				switch (intersection(*other.start, *other.end, *current.start, *current.end, first, second)) {
				case 2: intersections.push_back(second); intersections.push_back(first); break;
				case 1: if (last != 1 || first.z != 0) { intersections.push_back(first); } last = first.z; break;
				case 0: break;
				}
			}
			if (intersections.size() > 1) {
				std::sort(intersections.begin(), intersections.end(), [] (const Point &p, const Point &q) { return p.z < q.z; });
				size_t uniques = 1;
				for (size_t j = 1; j < intersections.size(); ++j) {
					if (
							std::fabs(intersections[j].x - intersections[j - 1].x) > 1e-12 ||
							std::fabs(intersections[j].y - intersections[j - 1].y) > 1e-12) {
						intersections[uniques++] = intersections[j];
					}
				}
				intersections.resize(uniques);
			}

			if (intersections.size() > 1) {
				if (intersections.size() % 2 == 1) {
					eslog::error("ESPRESO internal error: odd intersections.\n");
				}
				for (size_t j = 0; j < intersections.size(); j += 2) {
					lines.push_back({ intersections[j], intersections[j + 1] });
				}
			}
		}
	};
	std::vector<__line__> lines;
	addInnerLines(__line_iterator__(p), __line_iterator__(q), lines);
	addInnerLines(__line_iterator__(q), __line_iterator__(p), lines);

	auto equal = [] (Point &p, Point &q) {
		return std::fabs(p.x - q.x) < 1e-12 && std::fabs(p.y - q.y) < 1e-12;
	};

	while (lines.size()) {
		res.push_back({});
		res.back().push_back(lines.front().start);
		res.back().push_back(lines.front().end);
		lines.erase(lines.begin());
		size_t counter = 0;
		for (auto line = lines.begin(); counter < lines.size();) {
			if (line == lines.end()) {
				line = lines.begin();
			}
			if (equal(res.back().back(), line->start)) {
				res.back().push_back(line->end);
				line = lines.erase(line);
				counter = 0;
			} else {
				++line;
				++counter;
			}
			if (equal(res.back().back(), res.back().front())) {
				break;
			}
		}
		if (res.back().size() < 4 || !equal(res.back().back(), res.back().front())) {
			res.pop_back();
		} else {
			res.back().pop_back();
		}
	}
}

#define PERP_EPS 0.000001
#define dot2(a,b) ((a).x*(b).x + (a).y*(b).y)
#define dot3(a,b) ((a).x*(b).x + (a).y*(b).y + (a).z*(b).z)
#define perp(a,b) ((a).x*(b).y - (a).y*(b).x)

#define POLY_MAX 20 // maximum points in polygon
#define WEILERATHERTON_RESULTADD_MINLENGTH_RATIO 0.001
#define WEILERATHERTON_ENDLINE_MINLENGTH_RATIO 1e-6

// tells if 2D point p lies left comparing to oriented line e0->e1
// returns
//   > 0  p left  of e0 e1
//   ==0  p       on e0 e1
//   < 0  p right of e0 e1
// based on 3D vector product [(e1-e0).x, (e1-e0).y, 0] x [(p-e0).x, (p-e0).y, 0] = [0, 0, return value]
inline LeftOnRight isPointLeftToLine(Point p, Point e0, Point e1) {
	double where = perp(e1-e0,p-e0)/(e1-e0).length();  //(e1.x - e0.x) * (p.y - e0.y) - (p.x - e0.x) * (e1.y - e0.y);
	if (where > 0+PERP_EPS) {        return Left;
	} else if (where < 0-PERP_EPS) { return Right;
	} else {                         return On;
	}
}

bool is2dPointInLine(Point p, Point a0, Point a1) {
	Point a = a1 - a0;
	if (std::abs(a.x) >= std::abs(a.y)) {
		if (a0.x <= p.x && p.x <= a1.x) return true;
		if (a0.x >= p.x && p.x >= a1.x) return true;
	} else {
		if (a0.y <= p.y && p.y <= a1.y) return true;
		if (a0.y >= p.y && p.y >= a1.y) return true;
	}
	return false;
}

int intersect2DTwoLines(Point a0, Point a1, Point b0, Point b1, Point &i0, Point & i1) {

	double alen = (a1-a0).length();

	if        ((a0-b0).length() < alen*WEILERATHERTON_ENDLINE_MINLENGTH_RATIO){
		b0 = a0;
	} else if ((a1-b0).length() < alen*WEILERATHERTON_ENDLINE_MINLENGTH_RATIO){
		b0 = a1;
	}
	if        ((a0-b1).length() < alen*WEILERATHERTON_ENDLINE_MINLENGTH_RATIO){
		b1 = a0;
	} else if ((a1-b1).length() < alen*WEILERATHERTON_ENDLINE_MINLENGTH_RATIO){
		b1 = a1;
	}
	Point a    = a1 - a0;
	Point b    = b1 - b0;
	Point a0b0 = b0 - a0;
	double D    = perp(a,b);

	if (std::abs(D) < PERP_EPS) { // collinear lines
		if (std::abs(perp(a,a0b0)) > PERP_EPS || std::abs(perp(b,a0b0)) > PERP_EPS) { // two different parallel lines
			return 0;
		}
		double n2a = dot2(a,a);
		double n2b = dot2(b,b);
		if (n2a == 0 && n2b == 0) {
			if (a0 != b0) {
				return 0;
			}
			else {
				i0 = a0;
				return 1;
			}
		} else if (n2a == 0) {
			if (!is2dPointInLine(a0,b0,b1)) {
				return 0;
			} else {
				i0 = a0;
				return 1;
			}
		} else if (n2b == 0) {
			if (!is2dPointInLine(b0,a0,a1)) {
				return 0;
			} else {
				i0 = b0;
				return 1;
			}
		} else { // n2a != 0 && n2b != 0
			double t0, t1;
			Point a0b1 = b1 - a0;
			if (std::abs(a.x) >= std::abs(a.y)) {
				t0 = a0b0.x / a.x;
				t1 = a0b1.x / a.x;
			} else {
				t0 = a0b0.y / a.y;
				t1 = a0b1.y / a.y;
			}
			if (t0 > t1){
				double tmp = t0;
				t0 = t1;
				t1 = tmp;
			}
			if (t0 > 1 + WEILERATHERTON_ENDLINE_MINLENGTH_RATIO*alen || t1 < 0 - WEILERATHERTON_ENDLINE_MINLENGTH_RATIO*alen) {
				return 0;
			} else {
				t0 = (t0 < 0)? 0 : t0;
				t1 = (t1 > 1)? 1 : t1;
				if (t0 == t1){
					i0 = a0 + a*t0;
					return 1;
				} else {
					i0 = a0 + a*t0;
					i1 = a0 + a*t1;
					return 2;
				}
			}
		}
	} else { // not collinear lines
		float ta = -perp(b,a0b0) / D;
		if (ta < 0 - WEILERATHERTON_ENDLINE_MINLENGTH_RATIO*alen || ta > 1 + WEILERATHERTON_ENDLINE_MINLENGTH_RATIO*alen){
			return 0;
		}
		float tb = -perp(a,a0b0) / D;
		if (tb < 0 - WEILERATHERTON_ENDLINE_MINLENGTH_RATIO*alen || tb > 1 + WEILERATHERTON_ENDLINE_MINLENGTH_RATIO*alen){
			return 0;
		}
		i0 = a0 + a*ta;
		return 1;
	}
	return 0;
}

double getParamOfPointOnLine(Point p0, Point p1, Point p) {
	Point u = p1-p0;
	if (std::abs(u.x) < std::abs(u.y)) {
		return (p.y-p0.y)/u.y;
	} else {
		return (p.x-p0.x)/u.x;
	}
}

void clip(std::vector<Point> &p, std::vector<Point> &q, std::vector<std::vector<Point> > &res)
{
//	std::cout << std::endl << "-------------------------------------------------------------------------------------------------------------" << std::endl;
//	std::cout << "p:"; for (unsigned int i = 0; i < p.size(); i++) { std::cout << "(" << std::setprecision(18) << p[i].x << "," << std::setprecision(18) << p[i].y << ")-->";} std::cout << std::endl;
//	std::cout << "q:"; for (unsigned int i = 0; i < q.size(); i++) { std::cout << "(" << std::setprecision(18) << q[i].x << "," << std::setprecision(18) << q[i].y << ")-->";} std::cout << std::endl;
	// firstly try bounding boxes
	std::vector<Point> pBounds = {Point(DBL_MAX,DBL_MAX,0), Point(DBL_MIN,DBL_MIN,0)};
	std::vector<Point> qBounds = {Point(DBL_MAX,DBL_MAX,0), Point(DBL_MIN,DBL_MIN,0)};
	for (auto it = p.begin(); it != p.end(); it++){
		if (it->x < pBounds[0].x) pBounds[0].x = it->x;
		if (it->y < pBounds[0].y) pBounds[0].y = it->y;
		if (it->x > pBounds[1].x) pBounds[1].x = it->x;
		if (it->y > pBounds[1].y) pBounds[1].y = it->y;
	}
	for (auto it = q.begin(); it != q.end(); it++){
		if (it->x < qBounds[0].x) qBounds[0].x = it->x;
		if (it->y < qBounds[0].y) qBounds[0].y = it->y;
		if (it->x > qBounds[1].x) qBounds[1].x = it->x;
		if (it->y > qBounds[1].y) qBounds[1].y = it->y;
	}
	if (pBounds[0].x > qBounds[1].x || pBounds[1].x < qBounds[0].x || pBounds[0].y > qBounds[1].y || pBounds[1].y < qBounds[0].y) return;
	// then WeilerAtherton       https://www.youtube.com/watch?v=LCMyWFxeuro
	// should work for nonconvex - nonconvex clipping

	double typicalLength = std::min((p[0]-p[1]).length(),(q[0]-q[1]).length());
	std::vector<polyItem> pList(p.size());   for (auto it = pList.begin(); it != pList.end(); it++) { it->wn = 0;}
	std::vector<polyItem> qList(q.size());   for (auto it = qList.begin(); it != qList.end(); it++) { it->wn = 0;}
	std::vector<crossItem> cList;  cList.reserve(3*POLY_MAX);
//	std::cout << "WeilerAtherton Prepare  begins  -----------------------------------------------------------------------------------" << std::endl;
	for (unsigned int i = 0; i < p.size(); i++) {
		int j = (i+1) % p.size();
		// for all lines p[i]
		//   * compute all intersections of the line p[i]-->p[j] with all lines from q
		//   * decide if p[i] is inside q (winding number)
		Point xi = p[i];
		Point xj = p[j];
		for (unsigned int m = 0; m < q.size(); m++) {
			int n = (m+1) % q.size();
			Point xm = q[m];
			Point xn = q[n];
			// winding number of xi ------------------------------
			if (xm.y <= xi.y) { if (xn.y >  xi.y) { if (isPointLeftToLine(xi, xm, xn) == Left ) ++(pList[i].wn); } }
			else              { if (xn.y <= xi.y) { if (isPointLeftToLine(xi, xm, xn) == Right) --(pList[i].wn); } }
			// winding number of xk ------------------------------
			if (xi.y <= xm.y) { if (xj.y >  xm.y) { if (isPointLeftToLine(xm, xi, xj) == Left ) ++(qList[m].wn); } }
			else              { if (xj.y <= xm.y) { if (isPointLeftToLine(xm, xi, xj) == Right) --(qList[m].wn); } }
			// intersection xi-->xj with xk-->xl  ----------------
			Point intersection0, intersection1;
			double pt, qt;
			int nIntersects = intersect2DTwoLines(xi, xj, xm, xn, intersection0, intersection1);
			if (nIntersects > 0) {
//				std::cout << "  ";
//				std::cout << "p[" << std::setw(2) << i << "](" << std::setw(3) << p[i].x << "," << std::setw(3) << p[i].y << ")-->(" << std::setw(3) << p[j].x << "," << std::setw(3) << p[j].y << ") with ";
//				std::cout << "q[" << std::setw(2) << m << "](" << std::setw(3) << q[m].x << "," << std::setw(3) << q[m].y << ")-->(" << std::setw(3) << q[n].x << "," << std::setw(3) << q[n].y << ") ==> ";
//				std::cout << nIntersects << " intersections: (";
//				if (nIntersects < 2) { std::cout << "(" << std::setw(3) << intersection0.x << "," << std::setw(3) << intersection0.y << ")" << std::endl;}
//				else {                 std::cout << "(" << std::setw(3) << intersection0.x << "," << std::setw(3) << intersection0.y << "), "  << "(" << std::setw(3) << intersection1.x << "," << std::setw(3) << intersection1.y << ")" << std::endl;}
			}
			for (int inters = 0; inters < nIntersects; inters++) {
				if (inters > 0) {
					intersection0 = intersection1;
				}
				pt = getParamOfPointOnLine(xi, xj, intersection0);
				qt = getParamOfPointOnLine(xm, xn, intersection0);
				CrossTypes crt = decideCrossType(pt, qt, WEILERATHERTON_ENDLINE_MINLENGTH_RATIO);
				if (crt != None) {
					unsigned int i_(i), j_(j), k_((j+1) % p.size()), m_(m), n_(n), o_((m+2) % q.size());
					Point xi_(xi), xj_(xj), xk_(p[k_]), xm_(xm), xn_(xn), xo_(q[o_]);
					double pt_(pt), qt_(qt);
					switch (crt) {
					case  PDouble:
						if (pt < 0.5) {
							recomputeIndices(i_, p.size(), i_, j_, k_);  xk_ = xj_;  xj_ = xi_; xi_ = p[i_]; pt_ = 1.0;
						}
						break;
					case  QDouble:
						if (qt < 0.5) {
							recomputeIndices(m_, q.size(), m_, n_, o_);  xo_ = xn_;  xn_ = xm_; xm_ = q[m_]; qt_ = 1.0;
						}
						break;
					case PQDouble:
						if (pt < 0.5) {
							recomputeIndices(i_, p.size(), i_, j_, k_);  xk_ = xj_;  xj_ = xi_; xi_ = p[i_]; pt_ = 1.0;
						}
						if (qt < 0.5) {
							recomputeIndices(m_, q.size(), m_, n_, o_);  xo_ = xn_;  xn_ = xm_; xm_ = q[m_]; qt_ = 1.0;
						}
						break;
					default      : break;
					}
					// find if cross point exists in cList
					bool finded = false;
					unsigned int finded_ind(0);
					unsigned int tmp_cnt(0);
					for (auto it = cList.begin(); it != cList.end(); it++, tmp_cnt++) {
						if (it->i == i_ && it->m == m_ && std::abs(it->pPar-pt_) < 1e-4 && std::abs(it->qPar-qt_) < 1e-4) { finded = true; finded_ind = tmp_cnt;}
					}
					// end find
					NextLineTypes next;
					LeftOnRight xmLeftToxixj = isPointLeftToLine(xm_, xi_, xj_);
					LeftOnRight xnLeftToxixj = isPointLeftToLine(xn_, xi_, xj_);
					LeftOnRight xkLeftToxixj, xmLeftToxjxk, xnLeftToxjxk, xoLeftToxixj, xoLeftToxmxn, xoLeftToxjxk, xkLeftToxmxn, xkLeftToxnxo;
					switch (crt) {
					case   Single: // pt in (0,1)  and  qt in (0,1)   -- push_back (Single) Cross node
						next = decideSingleCrossNextLineType(xmLeftToxixj, xnLeftToxixj);
						if (finded) {
//							std::cout << " Node not added, because is already in cList" << std::endl;
						} else {
//							std::cout << " next=" << getNextLineTypeName(next) << std::endl;
							if (next != NLine) {  pushBackSingleNode(pList, qList, cList, intersection0, i , m , pt_, qt_, next); }
						}
						break;
					case  PDouble: // pt in {0,1}  and  qt in (0,1)   -- push_back  P Double Cross node
						xkLeftToxixj = isPointLeftToLine(xk_, xi_, xj_);
						xmLeftToxjxk = isPointLeftToLine(xm_, xj_, xk_);
						xnLeftToxjxk = isPointLeftToLine(xn_, xj_, xk_);
						next = decidePDoubleCrossNextLineType(xmLeftToxixj, xnLeftToxixj, xkLeftToxixj, xmLeftToxjxk, xnLeftToxjxk, dot2(xj_-xi_,xn-xm)<0);
						if (finded) {
							if (cList[finded_ind].type == Single) {
								alterToPDoubleNode( pList, qList, cList, finded_ind, intersection0, i_, m_, qt_, next);
							}else if (cList[finded_ind].type == QDouble) {
								alterToPQDoubleNode(pList, qList, cList, finded_ind, intersection0, i_, m_, next);
							} else {
//								std::cout << " Node not added, because is already in cList" << std::endl;
							}
						} else {
//							std::cout << " next=" << getNextLineTypeName(next) << std::endl;
							if (next != NLine) {  pushBackPDoubleNode(pList, qList, cList, intersection0, i_, m_, qt_, next); }
						}
						break;
					case  QDouble: // pt in (0,1)  and  qt in {0,1}   -- push_back  Q Double Cross node
						xoLeftToxixj = isPointLeftToLine(xo_, xi , xj );
						xoLeftToxmxn = isPointLeftToLine(xo_, xm_, xn_);
						next = decideQDoubleCrossNextLineType(xmLeftToxixj, xoLeftToxixj, xoLeftToxmxn, dot2(xj-xi,xn_-xm_)<0);
						if (finded) {
							if (cList[finded_ind].type == Single) {
								alterToQDoubleNode( pList, qList, cList, finded_ind, intersection0, i_, m_, pt_, next);
							}else if (cList[finded_ind].type == PDouble) {
								alterToPQDoubleNode(pList, qList, cList, finded_ind, intersection0, i_, m_, next);
							} else {
//								std::cout << " Node not added, because is already in cList" << std::endl;
							}
						} else {
//							std::cout << " next=" << getNextLineTypeName(next) << std::endl;
							if (next != NLine) {  pushBackQDoubleNode(pList, qList, cList, intersection0, i , m_, pt, next); }
						}
						break;
					case PQDouble: // pt in {0,1}  and  qt in {0,1}   -- push_back PQ Double Cross node
						xoLeftToxixj = isPointLeftToLine(xo_, xi_, xj_);
						xoLeftToxmxn = isPointLeftToLine(xo_, xm_, xn_);
						xkLeftToxixj = isPointLeftToLine(xk_, xi_, xj_);
						xoLeftToxjxk = isPointLeftToLine(xo_, xj_, xk_);
						xkLeftToxmxn = isPointLeftToLine(xk_, xm_, xn_);
						xkLeftToxnxo = isPointLeftToLine(xk_, xn_, xo_);
						xmLeftToxjxk = isPointLeftToLine(xm_, xj_, xk_);
						next = decidePQDoubleCrossNextLineType(xkLeftToxixj, xoLeftToxixj, xoLeftToxjxk, xoLeftToxmxn, xkLeftToxmxn, xkLeftToxnxo, xmLeftToxixj, xmLeftToxjxk, dot2(xj_-xi_,xn-xm)<0);
						if (finded) {
							if (cList[finded_ind].type == Single || cList[finded_ind].type == PDouble || cList[finded_ind].type == QDouble) {
								alterToPQDoubleNode(pList, qList, cList, finded_ind, intersection0, i_, m_, next);
							} else {
//								std::cout << " Node not added, because is already in cList" << std::endl;
							}
						} else {
//							std::cout << " next=" << getNextLineTypeName(next) << std::endl;
							if (next != NLine) {  pushBackPQDoubleNode(pList, qList, cList, intersection0, i_, m_, next); }
						}
						break;
					default: break;
					}
				}
			}
		}// end m
		// add xi to plist and sort incidence crosspoints list ----------------------------------
		sort(pList[i].inCross.begin(), pList[i].inCross.end(), [&](unsigned int  &a, unsigned int &b){ return (cList[a].pPar) < (cList[b].pPar);});
		if (!pList[i].insidesetted && pList[i].wn == 1) {
			pList[i].inside = true;
		}
	}// end i
	for (unsigned int i = 0; i < qList.size(); i++) {
		sort(qList[i].inCross.begin(), qList[i].inCross.end(), [&](unsigned int  &a, unsigned int &b){ return (cList[a].qPar) < (cList[b].qPar);});
		if (qList[i].wn == 1) {
			qList[i].inside = true;
		}
	}
	// cList ... fill pNext
	for (unsigned int i = 0; i < pList.size(); i++) {
		for (unsigned int j = 0; j < pList[i].inCross.size(); j++) {
			switch (cList[pList[i].inCross[j]].type) {
			case   Single:
			case  QDouble:
				if (j+1 < pList[i].inCross.size())  {
					cList[pList[i].inCross[j]].pNextType = CList;
					cList[pList[i].inCross[j]].pNext     = pList[i].inCross[j+1];
				} else {
					cList[pList[i].inCross[j]].pNextType = PList;
					cList[pList[i].inCross[j]].pNext     = (i+1) % p.size();
				}
				break;
			case  PDouble:
			case PQDouble:
				if (pList[(i+1) % pList.size()].inCross.size() > 0) {
					cList[pList[i].inCross[j]].pNextType = CList;
					cList[pList[i].inCross[j]].pNext     = pList[(i+1) % pList.size()].inCross[0];
				} else {
					cList[pList[i].inCross[j]].pNextType = PList;
					cList[pList[i].inCross[j]].pNext     = (i+2) % pList.size();
				}
				break;
			case     None:
				break;
			}
		}
	}
	// cList ... fill qNext
	for (unsigned int i = 0; i < qList.size(); i++) {
		for (unsigned int j = 0; j < qList[i].inCross.size(); j++) {
			switch (cList[qList[i].inCross[j]].type) {
			case   Single:
			case  PDouble:
				if (j+1 < qList[i].inCross.size()) {
					cList[qList[i].inCross[j]].qNext     = qList[i].inCross[j+1];
					cList[qList[i].inCross[j]].qNextType = CList;
				} else {
					cList[qList[i].inCross[j]].qNext     = (i+1) % qList.size();
					cList[qList[i].inCross[j]].qNextType = QList;
				}
				break;
			case  QDouble:
			case PQDouble:
				if (qList[(i+1) % qList.size()].inCross.size() > 0) {
					cList[qList[i].inCross[j]].qNextType = CList;
					cList[qList[i].inCross[j]].qNext     = qList[(i+1) % qList.size()].inCross[0];
				} else {
					cList[qList[i].inCross[j]].qNextType = QList;
					cList[qList[i].inCross[j]].qNext     = (i+2) % qList.size();
				}
				break;
			case     None:
				break;
			}
		}
	}
//	std::cout << "WeilerAtherton Prepare  ends  -----------------------------------------------------------------------------------" << std::endl;
	//-------------------------------------------------------------------------------------------
	// WEILER - ATHERTON
	// init
	WeilerAthertonState state(pList, qList, cList);
	int cnt = 0;
//	std::cout << "WeilerAtherton Loop begins -------------------------------------------------------------------------------------" << std::endl;
	if (cList.size() == 0) {
		unsigned int pInsides(0);
		unsigned int qInsides(0);
		for (unsigned int i = 0; i < pList.size(); i++) { if (pList[i].inside) pInsides++; }
		for (unsigned int i = 0; i < qList.size(); i++) { if (qList[i].inside) qInsides++; }
		if (pInsides == p.size()) {
			res.push_back(std::vector<Point>());
			for (auto it = p.begin(); it != p.end(); it++) {res.back().push_back(*it);}
			return;
		}
		if (qInsides == q.size()) {
			res.push_back(std::vector<Point>());
			for (auto it = q.begin(); it != q.end(); it++) {res.back().push_back(*it);}
			return;
		}
	} else {
		while(state.getNotEnd()) {
//			std::cout << "   processedNode: "   << state.getProcessedNode() << (state.isProcessingOutput()?"o":".");
//			std::cout << " outputStartNode: "   << state.getOutputStartNode() << " ";
//			std::cout << " pVisited: " << std::setw(2) << state.getPVisitedCountdown() << " ";
//			for (unsigned int i = 0; i < state.getPVisited().size(); i++) { std::cout << (state.getPVisited()[i]?"*":"."); }
//			std::cout << " cVisited: " << std::setw(2) << state.getCVisitedCountdown() << " ";
//			for (unsigned int i = 0; i < state.getCVisited().size(); i++) { std::cout << (state.getCVisited()[i]?"*":"."); }
			if (state.isProcessingOutput()) {       // ***** YES processing output **********************************************************
				if (state.wasVisited()) {           //       YES processing output  and YES visited *****************************************
					if (state.processedNodeIsOutputStartNode()) {
						state.endProcessOutput();
//						std::cout << "RESULT ";
//						for(unsigned int iii = 0; iii < res.back().size() ; iii++) {
//							std::cout << " (" << res.back()[iii].x << "," << res.back()[iii].y << ")-->";
//						}
						if (res.back().size() < 3) {
							res.pop_back();
//							std::cout << "ERASED" << std::endl;
						} else {
//							std::cout << "ADDED " << std::endl;
						}
						state.pNext();
					} else {
						state.setEnd();
//						std::cout << "Error:  isProcessingOutput() && wasVisited()" << std::endl;
					}
				} else {                            //       YES processing output  and NOT visited *****************************************
					switch (state.getProcessedNodeType()) {
					case PList:
						if (!state.isInside()) {
							state.setEnd();
//							std::cout << "Error:  isProcessingOutput() && !isInside() [PList]" << std::endl;
						} else {
							if (    ( res.back().back()-p[state.getProcessedNodeInd()]).length() > WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength) {
								if ((res.back().front()-p[state.getProcessedNodeInd()]).length() > WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength) {
//									std::cout << "     ||(" << ( res.back().back().x) << "," << ( res.back().back().y) << ")-(" << p[state.getProcessedNodeInd()].x << "," << p[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << ( res.back().back()-p[state.getProcessedNodeInd()]).length() << " > " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> push_back()   ";
									res.back().push_back(p[state.getProcessedNodeInd()]);
								} else {
//									std::cout << "     ||(" << (res.back().front().x) << "," << (res.back().front().y) << ")-(" << p[state.getProcessedNodeInd()].x << "," << p[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << (res.back().front()-p[state.getProcessedNodeInd()]).length() << " < " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> NOTHING       ";
								}
							} else {
//								std::cout <<     "     ||(" << ( res.back().back().x) << "," << ( res.back().back().y) << ")-(" << p[state.getProcessedNodeInd()].x << "," << p[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << ( res.back().back()-p[state.getProcessedNodeInd()]).length() << " < " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> NOTHING       ";
							}
							state.pNext();
						}
						break;
					case QList:
						if (!state.isInside()) {
							state.setEnd();
//							std::cout << "Error:  isProcessingOutput() && !isInside() [QList]" << std::endl;
						} else {
							if (    ( res.back().back()-q[state.getProcessedNodeInd()]).length() > WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength) {
								if ((res.back().front()-q[state.getProcessedNodeInd()]).length() > WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength) {
//									std::cout << "     ||(" << ( res.back().back().x) << "," << ( res.back().back().y) << ")-(" << q[state.getProcessedNodeInd()].x << "," << q[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << ( res.back().back()-q[state.getProcessedNodeInd()]).length() << " > " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> push_back()   ";
									res.back().push_back(q[state.getProcessedNodeInd()]);
								} else {
//									std::cout << "     ||(" << (res.back().front().x) << "," << (res.back().front().y) << ")-(" << q[state.getProcessedNodeInd()].x << "," << q[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << (res.back().front()-q[state.getProcessedNodeInd()]).length() << " < " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> NOTHING       ";
								}
							} else{
//								std::cout <<     "     ||(" << ( res.back().back().x) << "," << ( res.back().back().y) << ")-(" << q[state.getProcessedNodeInd()].x << "," << q[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << ( res.back().back()-q[state.getProcessedNodeInd()]).length() << " < " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> NOTHING       ";
							}
							state.qNext();
						}
						break;
					default :// CList
						if (    ( res.back().back()-Point(cList[state.getProcessedNodeInd()].x,cList[state.getProcessedNodeInd()].y,0)).length() > WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength) {
							if ((res.back().front()-Point(cList[state.getProcessedNodeInd()].x,cList[state.getProcessedNodeInd()].y,0)).length() > WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength) {
//								std::cout << "     ||(" << ( res.back().back().x) << "," << ( res.back().back().y) << ")-(" << cList[state.getProcessedNodeInd()].x << "," << cList[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << ( res.back().back()-MCVec2(cList[state.getProcessedNodeInd()].x,cList[state.getProcessedNodeInd()].y)).length() << " > " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> push_back()   ";
								res.back().push_back(Point(cList[state.getProcessedNodeInd()].x,cList[state.getProcessedNodeInd()].y,0));
							} else {
//								std::cout << "     ||(" << (res.back().front().x) << "," << (res.back().front().y) << ")-(" << cList[state.getProcessedNodeInd()].x << "," << cList[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << (res.back().front()-MCVec2(cList[state.getProcessedNodeInd()].x,cList[state.getProcessedNodeInd()].y)).length() << " < " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> NOTHING       ";
							}
						} else {
//							std::cout <<     "     ||(" << ( res.back().back().x) << "," << ( res.back().back().y) << ")-(" << cList[state.getProcessedNodeInd()].x << "," << cList[state.getProcessedNodeInd()].y <<")||=" << std::setprecision(4) << ( res.back().back()-MCVec2(cList[state.getProcessedNodeInd()].x,cList[state.getProcessedNodeInd()].y)).length() << " < " << std::setprecision(4) << WEILERATHERTON_RESULTADD_MINLENGTH_RATIO*typicalLength << "  -> NOTHING       ";
						}
						state.cNext();
						break;
					}
				}
			} else {                                // ***** NOT processing output **********************************************************
				if (state.wasVisited()) {           //       NOT processing output  and YES visited *****************************************
					state.pNext();
				} else {                            //       NOT processing output  and NOT visited *****************************************
					switch (state.getProcessedNodeType()) {
					case PList:
						if (state.isInside()) {
							state.startProcessOutput();
							res.push_back(std::vector<Point>());
							res.back().reserve(POLY_MAX);
							res.back().push_back(p[state.getProcessedNodeInd()]);
						}
						state.pNext();
						break;
					case QList:
						state.setEnd();
//						std::cout << "Error:  isfalse(state.processingOutput) && state.processedNode.inNodeType == QList" << std::endl;
						break;
					default :// CList
						if (state.canStartProcessOutput()) {
							state.startProcessOutput();
							res.push_back(std::vector<Point>());
							res.back().reserve(POLY_MAX);
							res.back().push_back(Point(cList[state.getProcessedNodeInd()].x,cList[state.getProcessedNodeInd()].y,0));
							state.cNext();
						} else {
							state.pNext();
						}
						break;
					}
				}
			}
			cnt++;
//			std::cout << cnt << std::endl;
			if (cnt > 80) state.setEnd();
		}
		if (state.isProcessingOutput()) {
			if (!res.empty() && res.back().size() < 3) { res.pop_back(); }
		}
	}
//	std::cout << "WeilerAtherton Loop  ends  -------------------------------------------------------------------------------------" << std::endl;


	//-------------------------------------------------------------------------------------------
	// OUTPUTS
//	std::cout << "plist----------------------" << std::endl;
//	for (unsigned int i = 0; i < p.size(); i++) {
//		//unsigned int j = (i+1) % p.size();
//		std::cout << "   " << std::setw(2) << i << pList[i] << std::endl;
//	}
//	std::cout << "qlist----------------------" << std::endl;
//	for (unsigned int i = 0; i < q.size(); i++) {
//		//unsigned int j = (i+1) % q.size();
//		std::cout << "   " << std::setw(2) << i << qList[i] << std::endl;
//	}
//	std::cout << "cross points--------------" << std::endl;
//	for (unsigned int i = 0; i < cList.size(); i++) {
//		std::cout << "   " << std::setw(2) << i << cList[i] << std::endl;
//	}
//	std::cout << "result--------------------" << std::endl;
//	for (unsigned int i = 0; i < res.size(); i++) {
//		std::cout << "   " << std::setw(2) << i << " ";
//		for (unsigned int j = 0; j < res[i].size(); j++) {
//			std::cout << "(" << std::setw(8) << res[i][j].x << "," << std::setw(8) << res[i][j].y << ")-->";
//		}
//		std::cout << std::endl;
//	}
}

}
}
