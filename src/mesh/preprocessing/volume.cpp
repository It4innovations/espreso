
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

#include <vector>

#include <fstream>
#include <sstream>
#include <iomanip>

#include "esinfo/meshinfo.h"

#include <unistd.h>

namespace espreso {
namespace mesh {

bool triangle_ray_intersect(Point p0, Point p1, Point v0, Point v1, Point v2);
bool edge_ray_intersect(Point p, Point v0, Point v1);

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

void store(Point voxels, const std::vector<int> &grid)
{
	std::ofstream casefile("volume.case");
	casefile
	<< "\n"
	<< "FORMAT\n"
	<< "type: ensight gold\n"
	<< "\n"
	<< "GEOMETRY\n"
	<< "model: volume.geo\n"
	<< "\n"
	<< "VARIABLE\n"
	<< "scalar per node: 1 VOLUME VOLUME.****\n"
	<< "\n"
	<< "TIME\n"
	<< "time set:               1\n"
	<< "number of steps:        1\n"
	<< "filename start numbers: 0\n"
	<< "filename increment:     1\n"
	<< "time values: 0\n";

	std::ofstream geo("volume.geo");
	geo
	<< "Output of our simple app for testing different ways of parallelization.\n"
	<< "You can open this file by majority of visualization tools (e.g., by ParaView: 'paraview heatflow.case')\n"
	<< "node id off\n"
	<< "element id off\n"
	<< "part\n"
	<< "         1\n"
	<< "2D uniform-elements (description line for part 1)\n"
	<< "block uniform\n";
	geo.width(10); geo << voxels.x;
	geo.width(10); geo << voxels.y;
	geo.width(10); geo << 1;
	geo << "\n";

	geo.precision(5);
	geo.setf(std::ios::scientific);
	geo.setf(std::ios::showpos);
	geo << 0. << "\n" << 0. << "\n" << 0. << "\n" << 1. / voxels.x << "\n" << -1. / voxels.y << "\n" << 0. << "\n";

	std::stringstream ss; ss << std::setw(4) << std::setfill('0') << 0;
	std::ofstream volume("VOLUME." + ss.str());

	volume
	<< "VOLUME\n"
	<< "part\n"
	<< "         1\n"
	<< "coordinates\n";

	volume.precision(5);
	volume.setf(std::ios::scientific);
	volume.setf(std::ios::showpos);

	for (esint r = 0; r < voxels.y; ++r) {
		for (esint c = 0; c < voxels.x; ++c) {
			volume << grid[r * voxels.x + c] << "\n";
		}
	}
}

void store3D(Point voxels, const std::vector<int> &grid, Point origin, double grid_offset)
{
	std::ofstream casefile("volume.case");
	casefile
			<< "\n"
			<< "FORMAT\n"
			<< "type: ensight gold\n"
			<< "\n"
			<< "GEOMETRY\n"
			<< "model: volume.geo\n"
			<< "\n"
			<< "VARIABLE\n"
			<< "scalar per node: 1 VOLUME VOLUME.****\n"
			<< "\n"
			<< "TIME\n"
			<< "time set:               1\n"
			<< "number of steps:        1\n"
			<< "filename start numbers: 0\n"
			<< "filename increment:     1\n"
			<< "time values: 0\n";

	std::ofstream geo("volume.geo");
	geo
			<< "Output of our simple app for testing different ways of parallelization.\n"
			<< "You can open this file by majority of visualization tools (e.g., by ParaView: 'paraview heatflow.case')\n"
			<< "node id off\n"
			<< "element id off\n"
			<< "part\n"
			<< "         1\n"
			<< "3D uniform-elements (description line for part 1)\n"
			<< "block uniform\n";
	geo.width(10); geo << voxels.x;
	geo.width(10); geo << voxels.y;
	geo.width(10); geo << voxels.z;
	geo << "\n";

	geo.precision(5);
	geo.setf(std::ios::scientific);
	geo.setf(std::ios::showpos);
	geo << origin.x << "\n" << origin.y << "\n" << origin.z << "\n" << grid_offset << "\n" << -grid_offset << "\n" << grid_offset << "\n";

	std::stringstream ss; ss << std::setw(4) << std::setfill('0') << 0;
	std::ofstream volume("VOLUME." + ss.str());

	volume
			<< "VOLUME\n"
			<< "part\n"
			<< "         1\n"
			<< "coordinates\n";

	volume.precision(5);
	volume.setf(std::ios::scientific);
	volume.setf(std::ios::showpos);

	for (esint t = 0; t < voxels.z; ++t) {
		for (esint r = 0; r < voxels.y; ++r) {
			for (esint c = 0; c < voxels.x; ++c) {
				volume << grid[t * voxels.y * voxels.x + r * voxels.x + c] << "\n";
			}
		}
	}
}

void computeVolumeIndices(ElementStore *elements, const NodeStore *nodes)
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
	std::vector<std::vector<_Point<int> > > vdata(info::env::OMP_NUM_THREADS);
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
								vdata[t].push_back(_Point<int>(x, y, z));
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
	std::vector<esint> distribution;
	std::vector<_Point<int> > data;
	utils::threadDistributionToFullDistribution(vdistribution);

	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++){
		distribution.reserve(distribution.size() + vdistribution[t].size());
		data.reserve(data.size() + vdata[t].size());

		distribution.insert(distribution.end(), vdistribution[t].begin(), vdistribution[t].end());
		data.insert(data.end(), vdata[t].begin(), vdata[t].end());
	}

	if(dim == 3){
		store3D(grid_size, grid, grid_start, grid_offset);
	} else {
		store(grid_size, grid);
	}

	elements->volumeIndices = new serializededata<esint, _Point<int> >(distribution, data);
	profiler::syncend("compute_volume_indices");
	eslog::checkpointln("MESH: VOLUME INDICES COMPUTED");
}

void OMcomputeVolumeIndices(ElementStore *elements, const NodeStore *nodes)
{
	profiler::syncstart("compute_volume_indices");

	double start = eslog::time();

	Point size(nodes->uniqInfo.max - nodes->uniqInfo.min), origin(nodes->uniqInfo.min);
	double step = (1.001 * std::max(std::max(size.x, size.y), size.z)) / info::ecf->output.volume_density;
	_Point<int> grid(std::floor(size.x / step) + 1, std::floor(size.y / step) + 1, std::floor(size.z / step) + 1);
	eslog::info(" == VOLUMETRIC GRID %55d x %5d x %5d == \n", grid.x, grid.y, grid.z);

	// set voxels to grid centers
	origin -= (Point(grid) * step - size) / 2;
	size = Point(grid) * step;

//	printf("[%f %f %f] -> [%f %f %f]\n", origin.x, origin.y, origin.z, origin.x + size.x, origin.y + size.y, origin.z + size.z);

	std::vector<std::vector<esint> > vdistribution(info::env::OMP_NUM_THREADS);
	std::vector<std::vector<_Point<int> > > vdata(info::env::OMP_NUM_THREADS);
	vdistribution[0].push_back(0);

	Point *coordinates = nodes->coordinates->datatarray().data();

	printf("init: %e\n", eslog::time() - start);
	start = eslog::time();

	double box = 0, angle = 0;
	size_t count = 0;

	if(info::mesh->dimension == 3) {
		#pragma omp parallel for
		for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
			auto epointers = elements->epointers->datatarray().begin(t);
			size_t eindex = 1;
			for (auto e = elements->nodes->cbegin(t); e != elements->nodes->cend(t); ++e, ++epointers, ++eindex) {

				double ss = eslog::time();
				// find BB of the element
				Point cmin = coordinates[e->front()], cmax = cmin;
				if (eindex == 0) printf("coo [%f %f %f]\n", coordinates[e->front()].x, coordinates[e->front()].y, coordinates[e->front()].z);
				for (auto n = e->begin() + 1; n != e->end(); ++n) {
					coordinates[*n].minmax(cmin, cmax);
					if (eindex == 0) printf("coo [%f %f %f]\n", coordinates[*n].x, coordinates[*n].y, coordinates[*n].z);
				}
				_Point<float> pmin = (cmin - origin) / size, pmax = (cmax - origin) / size;
				_Point<int> bmin(std::floor(pmin.x * grid.x), std::floor(pmin.y * grid.y), std::floor(pmin.z * grid.z));
				_Point<int> bmax(std::ceil(pmax.x * grid.x), std::ceil(pmax.y * grid.y), std::ceil(pmax.z * grid.z));
				if (eindex == 0) printf("[%f %f %f] -> [%f %f %f]\n", cmin.x, cmin.y, cmin.z, cmax.x, cmax.y, cmax.z);
				if (eindex == 0) printf("[%d %d %d] -> [%d %d %d]\n", bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z);
				box += eslog::time() - ss;
				ss = eslog::time();

				// loop through grid part points
				for (int x = bmin.x; x < bmax.x; ++x) {
					for (int y = bmin.y; y < bmax.y; ++y) {
						for (int z = bmin.z; z < bmax.z; ++z) {
							++count;
							const _Point<float> p = origin + _Point<float>(x, y, z) * step + .5 * step;
							if (eindex == 0) printf("[%f %f %f]\n", p.x, p.y, p.z);

							// loop through triangles
							float total_angle = 0;
							for (auto triangle = (*epointers)->triangles->begin(); triangle != (*epointers)->triangles->end(); ++triangle) {
								const _Point<float> v0 = coordinates[e->at(triangle->at(0))], v1 = coordinates[e->at(triangle->at(1))], v2 = coordinates[e->at(triangle->at(2))];
								if (v0 == p || v1 == p || v2 == p) {
									total_angle = 2 * M_PI;
									break;
								}

								total_angle += face_solid_angle_contribution(p, v0, v1, v2);
							}

							if (eindex == 0) printf("angle: %f\n", total_angle);
							// save element index if point is in the element
							if(fabs(total_angle) > 6.0) { // 2pi or 4pi -> inside
								vdata[t].push_back(_Point<int>(x, y, z));
							}
						}
					}
				}
				angle += eslog::time() - ss;
				vdistribution[t].push_back(vdata[t].size());
			}
		}
	} else {
		// TODO
	}
	printf("process: %e, box: %e, angle: %e, count: %lu, apv: %e\n", eslog::time() - start, box, angle, count, angle / count);
	start = eslog::time();

	utils::threadDistributionToFullDistribution(vdistribution);
	printf("points inside: %d\n", vdistribution.back().back());

	elements->volumeGrid = grid;
	elements->volumeOrigin = origin;
	elements->volumeSize = size;
	elements->volumeIndices = new serializededata<esint, _Point<int> >(vdistribution, vdata);

	printf("finish: %e\n", eslog::time() - start);

	profiler::syncend("compute_volume_indices");
	eslog::checkpointln("MESH: VOLUME INDICES COMPUTED");
}

bool triangle_ray_intersect(Point p0, Point p1, Point v0, Point v1, Point v2){
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

bool edge_ray_intersect(Point p, Point v0, Point v1){
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

}
}
