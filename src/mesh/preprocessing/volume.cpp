
#include "meshpreprocessing.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "wrappers/mpi/communication.h"

#include <vector>

#include <fstream>
#include <sstream>
#include <iomanip>

#include "esinfo/meshinfo.h"

#include <unistd.h>

namespace espreso {
namespace mesh {

static bool triangle_ray_intersect(Point p0, Point p1, Point v0, Point v1, Point v2){
	// triangle edge vectors and plane normal
	Point u = v1 - v0;
	Point v = v2 - v0;
	Point n = Point::cross(u, v);

	Point ray_dir = p1 - p0;
	Point w0 = p0 - v0;
	double a = -(n * w0);
	double b = n * ray_dir;

	double SMALL_NUM = 0.00000001; // ?
	if(fabs(b) < SMALL_NUM){ // ray is parallel to a plane
		if(a == 0){ // ray lies in triangle plane
			return false; //test if ray intersects triangle in 2D?
		} else { // ray is disjoint from triangle plane
			return false;
		}
	}

	// get intersect point of ray and triangle plane
	double r = a/b;
	if(r < 0.0){ // ray goes away from triangle
		return false;
	}

	Point i = p0 + ray_dir * r; // intersect point of ray and plane

	// is intersect point inside triangle
	double uu, uv, vv, wu, wv, d;
	Point w;
	uu = u*u;
	uv = u*v;
	vv = v*v;
	w = i - v0;
	wu = w*u;
	wv = w*v;
	d = uv * uv - uu*vv;

	// get and test parametric coords
	double s, t;
	s = (uv * wv - vv * wu) / d;
	if(s < 0.0 || s > 1.0){ // i is outside the triangle
		return false;
	}
	t = (uv * wu - uu * wv) / d;
	if(t < 0.0 || (s + t) > 1.0){ // i is outside the triangle
		return false;
	}

	return true; // i is in the triangle
}

static bool edge_ray_intersect(Point p, Point v0, Point v1){
	if(((v0.y <= p.y) && (v1.y > p.y)) // upward crossing
	   || ((v0.y > p.y) && (v1.y <= p.y))) { //downward crossing
		// edge-ray intersect
		float vt = (float)(p.y - v0.y) / (v1.y - v0.y);
		float x_intersect = v0.x + vt*(v1.x - v0.x);
		if(p.x < x_intersect){
			return true; //valid crossing right of p.x
		}
	}
	return false;
}

static inline float solid_angle(_Point<float> &a, _Point<float> &b, _Point<float> &c){
	a.normalize();
	b.normalize();
	c.normalize();

	float numer = _Point<float>::cross(a, b) * c;
	float denom = 1 + (a * b) + (b * c) + (c * a);

	float angle = 2 * atan2f(numer, denom);
	return fabsf(angle);
}

static inline float face_solid_angle_contribution(const _Point<float> &p, const _Point<float> &v0, const _Point<float> &v1, const _Point<float> &v2){
	_Point<float> vec_a = v0 - p;
	_Point<float> vec_b = v1 - p;
	_Point<float> vec_c = v2 - p;

	float angle = solid_angle(vec_a, vec_b, vec_c);

	_Point<float> u = v1 - v0;
	_Point<float> v = v2 - v0;
	_Point<float> n = _Point<float>::cross(u, v);

	_Point<float> p_vec = p - v0;
	float dot = n * p_vec;

	return dot > 0 ? angle : -angle;
}

static int ecode(const Element::CODE &code)
{
	switch (code) {
	case Element::CODE::POINT1:
		return 1;
	case Element::CODE::LINE2:
		return 3;
	case Element::CODE::LINE3:
		return 21;
	case Element::CODE::SQUARE4:
		return 9;
	case Element::CODE::SQUARE8:
		return 23;
	case Element::CODE::TRIANGLE3:
		return 5;
	case Element::CODE::TRIANGLE6:
		return 22;
	case Element::CODE::TETRA4:
		return 10;
	case Element::CODE::TETRA10:
		return 24;
	case Element::CODE::PYRAMID5:
		return 14;
	case Element::CODE::PYRAMID13:
		return 27;
	case Element::CODE::PRISMA6:
		return 13;
	case Element::CODE::PRISMA15:
		return 26;
	case Element::CODE::HEXA8:
		return 12;
	case Element::CODE::HEXA20:
		return 25;
	case Element::CODE::POLYHEDRON:
		return 42;
	default:
		return -1;
	}
}

void computeVolumeIndices(ElementStore *elements, const NodeStore *nodes)
{
	profiler::syncstart("compute_volume_indices");

	Point size(nodes->uniqInfo.max - nodes->uniqInfo.min), origin(nodes->uniqInfo.min);
	double step = (1.0001 * std::max(std::max(size.x, size.y), size.z)) / info::ecf->output.volume_density;
	_Point<short> grid(std::floor(size.x / step) + 1, std::floor(size.y / step) + 1, std::floor(size.z / step) + 1);
	eslog::info(" == GRID DENSITY %58d x %5d x %5d == \n", grid.x, grid.y, grid.z);
	Point epsilon(step / 100);

	// set voxels to grid centers
	origin -= (Point(grid) * step - size) / 2;
	origin += Point(step) / 2;
	size = Point(grid) * step;

	std::vector<std::vector<esint> > vdistribution(info::env::OMP_NUM_THREADS);
	std::vector<std::vector<_Point<short> > > vdata(info::env::OMP_NUM_THREADS);
	vdistribution[0].push_back(0);

	Point *coordinates = nodes->coordinates->datatarray().data();
	auto bbox = [&coordinates, &epsilon] (Element::CODE code, esint size, const esint *nodes, Point &min, Point &max) {
		min = max = coordinates[nodes[size - 1]]; // the last node is always coordinate
		PolyElement poly(Element::decode(code, size), nodes);
		for (esint n = 0; n < size; ++n) {
			if (poly.isNode(n)) {
				coordinates[nodes[n]].minmax(min, max);
			}
		}
		min -= epsilon;
		max += epsilon;
	};

	auto bboxIndices = [&origin, &size, &grid, &step] (const Point &min, const Point &max, _Point<short> &imin, _Point<short> &imax) {
		_Point<double> pmin = (min - origin) / size, pmax = (max - origin) / size;
		imin.x = std::floor(pmin.x * grid.x);
		imin.y = std::floor(pmin.y * grid.y);
		imin.z = std::floor(pmin.z * grid.z);
		imax.x = std::ceil(pmax.x * grid.x);
		imax.y = std::ceil(pmax.y * grid.y);
		imax.z = std::ceil(pmax.z * grid.z);

		_Point<double> p = origin + _Point<double>(imin) * step;
		int inside[3] = { 0, 0, 0 };
		for (int z = imin.z; z < imax.z; ++z, p.z += step) {
			if (min.z <= p.z && p.z <= max.z) {
				++inside[2];
			}
		}
		for (int y = imin.y; y < imax.y; ++y, p.y += step) {
			if (min.y <= p.y && p.y <= max.y) {
				++inside[1];
			}
		}
		for (int x = imin.x; x < imax.x; ++x, p.x += step) {
			if (min.x <= p.x && p.x <= max.x) {
				++inside[0];
			}
		}
		return inside[2] * inside[1] * inside[0];
	};

	auto cutConvex = [&] (const _Point<short> &imin, const _Point<short> &imax, const _Point<double> &v, const _Point<double> &n, std::vector<char> &isin) {
		_Point<double> p = origin + _Point<double>(imin) * step - v;
		int i = 0;
		for (short z = imin.z; z < imax.z; ++z, p.z += step) {
			for (short y = imin.y; y < imax.y; ++y, p.y += step) {
				for (short x = imin.x; x < imax.x; ++x, ++i, p.x += step) {
					isin[i / 8] &= ~((p * n < 0 ? 1 : 0) << (i % 8));
				}
				p.x = origin.x + imin.x * step - v.x;
			}
			p.y = origin.y + imin.y * step - v.y;
		}
	};

	auto cutConcave = [&] (const _Point<short> &imin, const _Point<short> &imax, const _Point<double> &v0, const _Point<double> &v1, const _Point<double> &v2, std::vector<double> &angle) {
		_Point<double> p = origin + _Point<double>(imin) * step;
		int i = 0;
		for (short z = imin.z; z < imax.z; ++z, p.z += step) {
			for (short y = imin.y; y < imax.y; ++y, p.y += step) {
				for (short x = imin.x; x < imax.x; ++x, ++i, p.x += step) {
					angle[i] += face_solid_angle_contribution(p, v0, v1, v2);
				}
				p.x = origin.x + imin.x * step;
			}
			p.y = origin.y + imin.y * step;
		}
	};

	if (info::mesh->dimension == 3) {
//		#pragma omp parallel for
		for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
			size_t eindex = 0;
			std::vector<_Point<short> > tdata; tdata.reserve(10 * elements->distribution.process.size);
			auto epointers = elements->epointers->datatarray().begin(t);
			for (auto e = elements->nodes->cbegin(t); e != elements->nodes->cend(t); ++e, ++epointers, ++eindex) {
				Point coomin, coomax;
				_Point<short> imin, imax;
				bbox((*epointers)->code, e->size(), e->data(), coomin, coomax);
				if (bboxIndices(coomin, coomax, imin, imax) == 0) {
					continue;
				}

				std::vector<char> isin((imax.x - imin.x) * (imax.y - imin.y) * (imax.z - imin.z) / 8 + 1, 255);
				std::vector<double> angle((imax.x - imin.x) * (imax.y - imin.y) * (imax.z - imin.z));
				switch ((*epointers)->code) {
				case Element::CODE::TETRA4: {
					// cut planes
					for (auto triangle = (*epointers)->triangles->begin(); triangle != (*epointers)->triangles->end(); ++triangle) {
						const _Point<double> &v = coordinates[e->at(triangle->at(0))];
						const _Point<double> n = _Point<double>::cross(coordinates[e->at(triangle->at(2))] - v, coordinates[e->at(triangle->at(1))] - v);
						cutConvex(imin, imax, v, n, isin);
					}
				} break;
				case Element::CODE::POLYHEDRON: {
					// compute the center of polygons since polygon can be non-planar!
					for (esint p = 0, pp = 1; p < e->front(); ++p, pp += e->at(pp) + 1) {
						Point pcenter;
						for (esint i = 1; i <= e->at(pp); ++i) {
							pcenter += coordinates[e->at(pp + i)];
						}
						pcenter /= e->at(pp);
						for (esint i = 0; i < e->at(pp); ++i) {
							const _Point<double> &v1 = coordinates[e->at(pp + 1 + i)];
							const _Point<double> &v2 = coordinates[e->at(pp + 1 + (i + 1) % e->at(pp))];
							cutConcave(imin, imax, pcenter, v1, v2, angle);
						}
					}
					for (size_t i = 0; i < angle.size(); ++i) {
						isin[i / 8] &= ~((fabs(angle[i]) < 6.0 ? 1 : 0) << (i % 8));
					}
				} break;
				default:
					for (auto face = (*epointers)->faces->begin(); face != (*epointers)->faces->end(); ++face) {
						// compute the center of each face since face can be non-planar!
						Point pcenter;
						for (size_t i = 0; i < face->size(); ++i) {
							pcenter += coordinates[e->at(face->at(i))];
						}
						pcenter /= face->size();
						for (size_t i = 0; i < face->size(); ++i) {
							const _Point<double> &v1 = coordinates[e->at(face->at(i))];
							const _Point<double> &v2 = coordinates[e->at(face->at((i + 1) % face->size()))];
							cutConcave(imin, imax, pcenter, v1, v2, angle);
						}
					}
					for (size_t i = 0; i < angle.size(); ++i) {
						isin[i / 8] &= ~((fabs(angle[i]) < 6.0 ? 1 : 0) << (i % 8));
					}
				}

				int i = 0;
				for (short z = imin.z; z < imax.z; ++z) {
					for (short y = imin.y; y < imax.y; ++y) {
						for (short x = imin.x; x < imax.x; ++x, ++i) {
							if (isin[i / 8] & (1 << (i % 8))) {
								tdata.push_back(_Point<short>(x, y, z));
							}
						}
					}
				}
				vdistribution[t].push_back(tdata.size());
			}
			vdata[t].swap(tdata);
		}
	} else {
		// TODO
	}

	utils::threadDistributionToFullDistribution(vdistribution);

	std::unique_lock<std::mutex> lk(info::mesh->voxelization.mutex);
	info::mesh->voxelization.cv.wait(lk, [] { return info::mesh->voxelization.counter == 0; });

	elements->volumeGrid = grid;
	elements->volumeOrigin = origin;
	elements->volumeSize = size;
	if (elements->volumeIndices) {
		delete elements->volumeIndices;
	}
	elements->volumeIndices = new serializededata<esint, _Point<short> >(vdistribution, vdata);

	profiler::syncend("compute_volume_indices");
	eslog::checkpointln("MESH: VOLUME INDICES COMPUTED");

	size_t totalHits = vdistribution.back().back();
	Communication::allReduce(&totalHits, nullptr, 1, MPITools::getType(totalHits).mpitype, MPI_SUM);
	eslog::info(" == NUMBER OF NON-EMPTY VOXELS %50lu [%6.2f%%] == \n", totalHits, 100.0 * totalHits / ((size_t)grid.x * grid.y * grid.z));

	eslog::checkpointln("MESH: VOLUME INDICES SYNCHRONIZATION");
}

void computeVolumeIndicesOM(ElementStore *elements, const NodeStore *nodes)
{
	profiler::syncstart("compute_volume_indices");

	Point size(nodes->uniqInfo.max - nodes->uniqInfo.min), origin(nodes->uniqInfo.min);
	double step = (1.0001 * std::max(std::max(size.x, size.y), size.z)) / info::ecf->output.volume_density;
	_Point<short> grid(std::floor(size.x / step) + 1, std::floor(size.y / step) + 1, std::floor(size.z / step) + 1);
	eslog::info(" == GRID DENSITY %58d x %5d x %5d == \n", grid.x, grid.y, grid.z);

	// set voxels to grid centers
	origin -= (Point(grid) * step - size) / 2;
	origin += Point(step) / 2;
	size = Point(grid) * step;

	std::vector<std::vector<esint> > vdistribution(info::env::OMP_NUM_THREADS);
	std::vector<std::vector<_Point<short> > > vdata(info::env::OMP_NUM_THREADS);
	vdistribution[0].push_back(0);

	Point *coordinates = nodes->coordinates->datatarray().data();

	if(info::mesh->dimension == 3) {
//		#pragma omp parallel for
		for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
			size_t eclose = 0;
			Point ppmax(10, 10, 10), ppboundmin, ppboundmax;
			size_t eindex = 0;
			std::vector<_Point<short> > tdata; tdata.reserve(10 * elements->distribution.process.size);
			auto epointers = elements->epointers->datatarray().begin(t);
			for (auto e = elements->nodes->cbegin(t); e != elements->nodes->cend(t); ++e, ++epointers, ++eindex) {
				int cc = 0;
				Point coomin = coordinates[e->back()], coomax = coomin, center;
				PolyElement poly(Element::decode((*epointers)->code, e->size()), e->begin());
				for (size_t n = 0; n < e->size(); ++n) {
					if (poly.isNode(n)) {
						center += coordinates[e->at(n)];
						coordinates[e->at(n)].minmax(coomin, coomax);
						++cc;
					}
				}
				center /= cc;
				Point dist = origin + _Point<double>(1, 0, 0) * step;
				if ((center - dist).length() < ppmax.length()) {
					ppmax = center;
					ppboundmin = coomin;
					ppboundmax = coomax;
					eclose = eindex;
				}
				if (coomin.x <= 0.999600 && 0.999600 <= coomax.x) {
					if (coomin.y <= -2.000200 && -2.000200 <= coomax.y) {
						if (coomin.z <= 1.999800 && 1.999800 <= coomax.z) {
							printf("adept %d %lu [%f %f %f] [%f %f %f]\n", info::mpi::rank, eindex, coomin.x, coomin.y, coomin.z, coomax.x, coomax.y, coomax.z);
						}
					}
				}

				_Point<double> pmin = (coomin - origin) / size, pmax = (coomax - origin) / size;
				_Point<short> bmin(std::floor(pmin.x * grid.x), std::floor(pmin.y * grid.y), std::floor(pmin.z * grid.z));
				_Point<short> bmax(std::ceil(pmax.x * grid.x), std::ceil(pmax.y * grid.y), std::ceil(pmax.z * grid.z));

				std::vector<char> isin((bmax.x - bmin.x) * (bmax.y - bmin.y) * (bmax.z - bmin.z) / 8 + 1, 255);
				auto cut = [&] (const _Point<double> &v, const _Point<double> &n) {
					_Point<double> p = origin + _Point<double>(bmin) * step - v;
					for (short z = bmin.z, i = 0; z < bmax.z; ++z, p.z += step) {
						for (short y = bmin.y; y < bmax.y; ++y, p.y += step) {
							for (short x = bmin.x; x < bmax.x; ++x, ++i, p.x += step) {
								isin[i / 8] &= ~((p * n < 0 ? 1 : 0) << (i % 8));
							}
							p.x = origin.x + bmin.x * step - v.x;
						}
						p.y = origin.y + bmin.y * step - v.y;
					}
				};

				if ((*epointers)->code == Element::CODE::POLYHEDRON) {
					bool out = true;
					_Point<double> p = origin + _Point<double>(bmin) * step;
					for (short z = bmin.z; z < bmax.z; ++z, p.z += step) {
						if (coomin.z <= p.z && p.z <= coomax.z) {
							out = false;
						}
					}
					for (short y = bmin.y; y < bmax.y; ++y, p.y += step) {
						if (coomin.y <= p.y && p.y <= coomax.y) {
							out = false;
						}
					}
					for (short x = bmin.x; x < bmax.x; ++x, p.x += step) {
						if (coomin.x <= p.x && p.x <= coomax.x) {
							out = false;
						}
					}
					if (out) {
						std::fill(isin.begin(), isin.end(), 0);
					}
					for (esint p = 0, pp = 1; p < e->front(); ++p, pp += e->at(pp) + 1) {
						// all polygons are in plane -> it is enough to check single triangle only
						// but, we need to check if the points are not at the line
						const _Point<double> &v = coordinates[e->at(pp + 1)];
						_Point<double> n2, n1 = _Point<double>::cross(coordinates[e->at(pp + 2)] - v, coordinates[e->at(pp + 3)] - v);
						if (e->at(pp) > 3) {
							n2 = _Point<double>::cross(coordinates[e->at(pp + 2)] - v, coordinates[e->at(pp + 4)] - v);
						}
						if (n1.x * n1.x + n1.y * n1.y + n1.z * n1.z < n2.x * n2.x + n2.y * n2.y + n2.z * n2.z) {
							cut(v, n2);
						} else {
							cut(v, n1);
						}
					}
				} else {
					for (auto triangle = (*epointers)->triangles->begin(); triangle != (*epointers)->triangles->end(); ++triangle) {
						const _Point<double> &v = coordinates[e->at(triangle->at(0))];
						const _Point<double> n = _Point<double>::cross(coordinates[e->at(triangle->at(2))] - v, coordinates[e->at(triangle->at(1))] - v);
						cut(v, n);
					}
				}
				for (short z = bmin.z, i = 0; z < bmax.z; ++z) {
					for (short y = bmin.y; y < bmax.y; ++y) {
						for (short x = bmin.x; x < bmax.x; ++x, ++i) {
							if (isin[i / 8] & (1 << (i % 8))) {
								printf("%d %d %d\n", x, y, z);
								tdata.push_back(_Point<short>(x, y, z));
							}
						}
					}
				}
//				if (info::mpi::rank == 2 && eindex == 3798) {
//					size_t voxx = tdata.size() - vdistribution[t].back();
//					std::ofstream os(std::to_string(eindex) + ".vtk");
//					os << "# vtk DataFile Version 2.0\n";
//					os << "element with voxels\n";
//					os << "ASCII\n";
//					os << "DATASET UNSTRUCTURED_GRID\n";
//					os << "POINTS " << e->size() + voxx << " float\n";
//					for (auto n = e->begin(); n != e->end(); ++n) {
//						os << coordinates[*n].x << " " << coordinates[*n].y << " " << coordinates[*n].z << "\n";
//					}
//					_Point<double> p = origin + _Point<double>(bmin) * step + .5 * step;
//					for (int z = bmin.z, i = 0; z < bmax.z; ++z, p.z += step) {
//						for (int y = bmin.y; y < bmax.y; ++y, p.y += step) {
//							for (int x = bmin.x; x < bmax.x; ++x, ++i, p.x += step) {
//								if (isin[i / 8] & (1 << (i % 8))) {
//									os << p.x << " " << p.y << " " << p.z << "\n";
//								}
//							}
//							p.x = origin.x + bmin.x * step + .5 * step;
//						}
//						p.y = origin.y + bmin.y * step + .5 * step;
//					}
//					os << "CELLS " << 1 + voxx << " " << e->size() + 1 + 2 * voxx << "\n";
//					os << e->size();
//					PolyElement poly(Element::decode((*epointers)->code, e->size()), e->begin());
//					for (size_t n = 0; n < e->size(); ++n) {
//						if (poly.isNode(n)) {
//							os << " " << n;
//						} else {
//							os << " " << e->at(n);
//						}
//					}
//					os << "\n";
//					for (size_t i = 0; i < voxx; ++i) {
//						os << "1 " << i + e->size() << "\n";
//					}
//					os << "CELL_TYPES " << voxx + 1 << "\n";
//					os << ecode((*epointers)->code) << "\n";
//					for (size_t i = 0; i < voxx; ++i) {
//						os << "1\n";
//					}
//					os.close();
//
//					printf(" << %lu [%d]\n", eindex, info::mpi::rank);
////					PolyElement polyx(Element::decode((*epointers)->code, e->size()), e->begin());
////					for (size_t n = 0; n < e->size(); ++n) {
////						if (polyx.isNode(n)) {
////							printf("%d ", nodes->IDs->datatarray()[e->at(n)]);
////						} else {
////							printf("%d ", e->at(n));
////						}
////					}
////					printf("\n");
//				}
				vdistribution[t].push_back(tdata.size());
			}
			vdata[t].swap(tdata);
			_Point<double> p100 = origin + _Point<double>(1, 0, 0) * step;
			printf("closest: %lu [%f %f %f] [%f %f %f] [%f %f %f]\n", eclose, p100.x, p100.y, p100.z, ppboundmin.x, ppboundmin.y, ppboundmin.z, ppboundmax.x, ppboundmax.y, ppboundmax.z);
		}
	} else {
		// TODO
	}

	utils::threadDistributionToFullDistribution(vdistribution);

	elements->volumeGrid = grid;
	elements->volumeOrigin = origin;
	elements->volumeSize = size;
	elements->volumeIndices = new serializededata<esint, _Point<short> >(vdistribution, vdata);

	profiler::syncend("compute_volume_indices");
	eslog::checkpointln("MESH: VOLUME INDICES COMPUTED");

	esint totalHits = vdistribution.back().back();
	Communication::allReduce(&totalHits, nullptr, 1, MPITools::getType(totalHits).mpitype, MPI_SUM);
	eslog::info(" == NUMBER OF NON-EMPTY VOXELS %50d [%5.2f%%] == \n", totalHits, 100.0 * totalHits / (grid.x * grid.y * grid.z));

	eslog::checkpointln("MESH: VOLUME INDICES SYNCHRONIZATION");
}

void computeVolumeIndicesMH(ElementStore *elements, const NodeStore *nodes)
{
	profiler::syncstart("compute_volume_indices");

	// find BB of the mesh
	const Point &p_min_max = nodes->coordinates->datatarray()[*(elements->nodes->cbegin()->begin())];
	Point mesh_min = Point(p_min_max.x, p_min_max.y, p_min_max.z);
	Point mesh_max = Point(p_min_max.x, p_min_max.y, p_min_max.z);
	for (auto e = elements->nodes->cbegin(); e != elements->nodes->cend(); ++e) {
		for (auto n = e->begin(); n != e->end(); ++n) {
			Point &p = nodes->coordinates->datatarray()[*n];
			p.minmax(mesh_min, mesh_max);
		}
	}

	Point mesh_min_global;
	Point mesh_max_global;
	MPI_Allreduce(&mesh_min, &mesh_min_global, 3, MPI_DOUBLE, MPI_MIN, info::mpi::comm);
	MPI_Allreduce(&mesh_max, &mesh_max_global, 3, MPI_DOUBLE, MPI_MAX, info::mpi::comm);

	// grid setting
	int grid_size_x = info::ecf->output.volume_density;
	double grid_offset = (mesh_max_global.x - mesh_min_global.x)/(double)(grid_size_x - 1);
	printf("offset: %f\n", grid_offset);
	Point grid_size = Point(grid_size_x, 
						(int)((mesh_max_global.y - mesh_min_global.y)/grid_offset + 2.0),
						(int)((mesh_max_global.z - mesh_min_global.z)/grid_offset + 2.0));
	int dim = info::mesh->dimension;
	std::vector<int> grid;
	int grid_init_value = -1;
	if(dim == 3){
		grid.resize(grid_size.x * grid_size.y * grid_size.z);
		fill(grid.begin(), grid.end(), grid_init_value);
	} else { // dim == 2
		grid.resize(grid_size.x * grid_size.y);
		fill(grid.begin(), grid.end(), grid_init_value);
	}

	Point grid_start = Point(mesh_min_global.x, mesh_max_global.y, mesh_min_global.z);

	std::vector<std::vector<esint> > vdistribution(info::env::OMP_NUM_THREADS);
	std::vector<std::vector<_Point<short> > > vdata(info::env::OMP_NUM_THREADS);
	vdistribution[0].push_back(0);

	// elements cycle
	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
	auto epointers = elements->epointers->cbegin(t);
		for (auto e = elements->nodes->cbegin(t); e != elements->nodes->cend(t); ++e, ++epointers) {
		
			// find BB of the element
			const Point &p_min_max = nodes->coordinates->datatarray()[*(e->begin())];
			Point el_min = Point(p_min_max.x, p_min_max.y, p_min_max.z);
			Point el_max = Point(p_min_max.x, p_min_max.y, p_min_max.z);

			for (auto n = e->begin(); n != e->end(); ++n) {
				Point &p = nodes->coordinates->datatarray()[*n];
				p.minmax(el_min, el_max);
			}

			// grid part in BB
			int min_x_inx = (el_min.x - grid_start.x)/grid_offset;
			int min_y_inx = (grid_start.y - el_max.y)/grid_offset;
			int min_z_inx = (el_min.z - grid_start.z)/grid_offset;
			int max_x_inx = (el_max.x - grid_start.x)/grid_offset;
			int max_y_inx = (grid_start.y - el_min.y)/grid_offset;
			int max_z_inx = (el_max.z - grid_start.z)/grid_offset;

			// loop through grid part points
			for(int x = min_x_inx; x <= max_x_inx; x++){
				for(int y = min_y_inx; y <= max_y_inx; y++){
					for(int z = min_z_inx; z <= max_z_inx; z++){
						// skip point, which is already assigned to some element
						int grid_inx_1d = z*grid_size.y*grid_size.x + y*grid_size.x + x;
						if(grid[grid_inx_1d] != -1){
							continue;
						}

						// test point in element
						Point p = Point(grid_start.x + x*grid_offset, grid_start.y - y*grid_offset, grid_start.z + z*grid_offset);

						if(dim == 3){
							// loop through triangles
							double total_angle = 0;
							bool is_p_triangle_vertex = false;

							auto triangles = epointers->at(0)->triangles;
							for (auto triangle = triangles->begin(); triangle != triangles->end(); ++triangle) {

								Point v0 = nodes->coordinates->datatarray()[e->at(*(triangle->begin()))];
								Point v1 = nodes->coordinates->datatarray()[e->at(*(triangle->begin() + 1))];
								Point v2 = nodes->coordinates->datatarray()[e->at(*(triangle->begin() + 2))];

								if(v0 == p || v1 == p || v2 == p){
									is_p_triangle_vertex = true;
									break;
								}

								total_angle += face_solid_angle_contribution(p, v0, v1, v2);
							}

							// save element index if point is in the element
							if(is_p_triangle_vertex || fabs(total_angle) > 6.0){ // 2pi or 4pi -> inside
								grid[grid_inx_1d] = 0;
								vdata[t].push_back(_Point<short>(x, y, z));
							}

						} else { // dim == 2
							int cn = 0; // cn

							// loop through edges
							for (auto n = e->begin(); n != e->end() - 1; ++n) {
								Point &v0 = nodes->coordinates->datatarray()[*n];
								Point &v1 = nodes->coordinates->datatarray()[*(n+1)];

								if(edge_ray_intersect(p, v0, v1)){
									cn++;
								}
							}

							// last edge
							Point &v0 = nodes->coordinates->datatarray()[*(e->end() - 1)];
							Point &v1 = nodes->coordinates->datatarray()[*(e->begin())];
							if(edge_ray_intersect(p, v0, v1)){
								cn++;
							}

							// save polygon index if point is in polygon
							if(cn%2){ // cn is odd
								grid[z*grid_size.y*grid_size.x + y*grid_size.x + x] = 0; 
							}
						}
					}
				}
			}
			vdistribution[t].push_back(vdata[t].size());
		}
	}
	// combine thread data
	utils::threadDistributionToFullDistribution(vdistribution);

//	elements->volumeGrid = grid;
//	elements->volumeOrigin = origin;
//	elements->volumeSize = size;
	elements->volumeIndices = new serializededata<esint, _Point<short> >(vdistribution, vdata);
	profiler::syncend("compute_volume_indices");
	eslog::checkpointln("MESH: VOLUME INDICES COMPUTED");
}

}
}
