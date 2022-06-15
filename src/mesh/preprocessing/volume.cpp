
#include "meshpreprocessing.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
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
double face_solid_angle_contribution(Point p, Point v0, Point v1, Point v2);
double solid_angle(Point a, Point b, Point c);

void store(const esint &voxels, const std::vector<int> &grid)
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
	geo.width(10); geo << voxels;
	geo.width(10); geo << voxels;
	geo.width(10); geo << 1;
	geo << "\n";

	geo.precision(5);
	geo.setf(std::ios::scientific);
	geo.setf(std::ios::showpos);
	geo << 0. << "\n" << 0. << "\n" << 0. << "\n" << 1. / voxels << "\n" << -1. / voxels << "\n" << 0. << "\n";

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

	for (esint r = 0; r < voxels; ++r) {
		for (esint c = 0; c < voxels; ++c) {
			volume << grid[r * voxels + c] << "\n";
		}
	}
}

void store3D(const esint &voxels, const std::vector<int> &grid, Point origin, Point grid_offset)
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
	geo.width(10); geo << voxels;
	geo.width(10); geo << voxels;
	geo.width(10); geo << voxels;
	geo << "\n";

	geo.precision(5);
	geo.setf(std::ios::scientific);
	geo.setf(std::ios::showpos);
	geo << origin.x << "\n" << origin.y << "\n" << origin.z << "\n" << grid_offset.x << "\n" << -grid_offset.y << "\n" << grid_offset.z << "\n";

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

	for (esint t = 0; t < voxels; ++t) {
		for (esint r = 0; r < voxels; ++r) {
			for (esint c = 0; c < voxels; ++c) {
				volume << grid[t * voxels * voxels + r * voxels + c] << "\n";
			}
		}
	}
}

void computeVolumeIndices(ElementStore *elements, const NodeStore *nodes)
{
	profiler::syncstart("compute_volume_indices");

	// uniform grid
	esint grid_size = info::ecf->output.volume_density;
	int voxels = grid_size - 1;
	int dim = info::mesh->dimension;
	std::vector<int> grid;
	int grid_init_value = -1;
	if(dim == 3){
		grid.resize(grid_size * grid_size * grid_size);
		fill(grid.begin(), grid.end(), grid_init_value);
	} else { // dim == 2
		grid.resize(grid_size * grid_size);
		fill(grid.begin(), grid.end(), grid_init_value);
	}

	// find BB of the mesh
	const Point &p_min_max = nodes->coordinates->datatarray()[*(elements->nodes->cbegin()->begin())];
	Point mesh_min = Point(p_min_max.x, p_min_max.y, p_min_max.z);
	Point mesh_max = Point(p_min_max.x, p_min_max.y, p_min_max.z);
	for (auto e = elements->nodes->cbegin(); e != elements->nodes->cend(); ++e) {
		for (auto n = e->begin(); n != e->end(); ++n) {
			Point &p = nodes->coordinates->datatarray()[*n];
			//printf("%f %f %f\n", p.x, p.y, p.z);

			p.minmax(mesh_min, mesh_max);
		}
	}
	printf("mesh min: %f %f %f\n", mesh_min.x, mesh_min.y, mesh_min.z);
	printf("mesh max: %f %f %f\n", mesh_max.x, mesh_max.y, mesh_max.z);

	// grid setting
	Point grid_start = Point(mesh_min.x, mesh_max.y, mesh_min.z);
	Point grid_end = Point(mesh_max.x, mesh_min.y, mesh_max.z);
	Point grid_offset = Point((grid_end.x - grid_start.x)/voxels,
							(grid_start.y - grid_end.y)/voxels, 
							(grid_end.z - grid_start.z)/voxels);
	printf("offset: %f %f %f\n", grid_offset.x, grid_offset.y, grid_offset.z);

	// debug
	// FILE * pFile;
   	// pFile = fopen ("volume_debug.txt","w");
	// fprintf(pFile, "x;y;z;x_inx;y_inx;z_inx;element_inx\n");

	// FILE * pFile_element;
   	// pFile_element = fopen ("./debug_files/point_element.txt","w");

	//FILE * pFile_hist;
   	//pFile_hist = fopen ("./debug_files/angle_results.txt","w");

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
			int min_x_inx = (el_min.x - grid_start.x)/grid_offset.x;
			int min_y_inx = (grid_start.y - el_max.y)/grid_offset.y;
			int min_z_inx = (el_min.z - grid_start.z)/grid_offset.z;
			int max_x_inx = (el_max.x - grid_start.x)/grid_offset.x;
			int max_y_inx = (grid_start.y - el_min.y)/grid_offset.y;
			int max_z_inx = (el_max.z - grid_start.z)/grid_offset.z;

			// loop through grid part points
			for(int x = min_x_inx; x <= max_x_inx; x++){
				for(int y = min_y_inx; y <= max_y_inx; y++){
					for(int z = min_z_inx; z <= max_z_inx; z++){
						// skip point, which is already assigned to some element
						int grid_inx_1d = z*grid_size*grid_size + y*grid_size + x;
						if(grid[grid_inx_1d] != -1){
							continue;
						}

						// test point in element
						Point p = Point(grid_start.x + x*grid_offset.x, grid_start.y - y*grid_offset.y, grid_start.z + z*grid_offset.z);

						if(dim == 3){
							// loop through faces
							//printf("faces:\n");
							//Point p1 = Point(el_max.x, p.y, p.z); // cn
							double total_angle = 0;
							bool is_p_face_vertex = false;

							auto faces = epointers->at(0)->faces;
							for (auto face = faces->begin(); face != faces->end(); ++face) {

								int num_verts = face->end() - face->begin();
								//printf("num face verts: %d\n", num_verts);

								if(num_verts == 3 || num_verts == 6){ // triangle
									Point v0 = nodes->coordinates->datatarray()[e->at(*(face->begin()))];
									Point v1 = nodes->coordinates->datatarray()[e->at(*(face->begin() + 1))];
									Point v2 = nodes->coordinates->datatarray()[e->at(*(face->begin() + 2))];

									if(v0 == p || v1 == p || v2 == p){
										is_p_face_vertex = true;
										break;
									}

									total_angle += face_solid_angle_contribution(p, v0, v1, v2);

								} else { // rectangle
									Point v0 = nodes->coordinates->datatarray()[e->at(*(face->begin()))];
									Point v1 = nodes->coordinates->datatarray()[e->at(*(face->begin() + 1))];
									Point v2 = nodes->coordinates->datatarray()[e->at(*(face->begin() + 2))];
									Point v3 = nodes->coordinates->datatarray()[e->at(*(face->begin() + 3))];
									
									if(v0 == p || v1 == p || v2 == p || v3 == p){
										is_p_face_vertex = true;
										break;
									}

									total_angle += face_solid_angle_contribution(p, v0, v1, v2);
									total_angle += face_solid_angle_contribution(p, v0, v2, v3);

									// debug
									//if(x==3 && y==87 && z==82){
										// fprintf(pFile_element, "%f;%f;%f\n", v0.x, v0.y, v0.z);
										// fprintf(pFile_element, "%f;%f;%f\n", v1.x, v1.y, v1.z);
										// fprintf(pFile_element, "%f;%f;%f\n", v2.x, v2.y, v2.z);
										// fprintf(pFile_element, "%f;%f;%f\n", v3.x, v3.y, v3.z);
										// printf("point %f %f %f\n", p.x, p.y, p.z);
									//}

								}
							}

							//debug
							//fprintf(pFile_hist, "%f\n", fabs(total_angle));

							// save element index if point is in the element
							if(is_p_face_vertex || fabs(total_angle) > 6.0){ // 2pi or 4pi -> inside
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
								grid[z*grid_size*grid_size + y*grid_size + x] = 0; 
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
	printf("here");


	//debug
	// for(int i = 0; i < grid.size(); i++){
	// 	int z = i/(grid_size*grid_size);
	// 	int y = (i-z*grid_size*grid_size)/grid_size;
	// 	int x = i - z*grid_size*grid_size - y*grid_size;	

	// 	fprintf(pFile, "%f;%f;%f;%d;%d;%d;%d\n", grid_start.x + x*grid_offset.x, grid_start.y - y*grid_offset.y, grid_start.z + z*grid_offset.z, x, y, z, grid[i]); // fprintf(pFile, "x;y;z;x_inx;y_inx;z_inx;element_inx;face_vertex;total_angle\n");
	// }	
	// fclose(pFile);
	// fclose(pFile_element);
	// fclose(pFile_hist);

	//grid[81*grid_size*grid_size + 87*grid_size + 3] = 0; // x==3 && y==87 && z==81

	if(dim == 3){
		store3D(grid_size, grid, grid_start, grid_offset);
	} else {
		store(grid_size, grid);
	}

	elements->volumeIndices = new serializededata<esint, _Point<int> >(distribution, data);
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

double face_solid_angle_contribution(Point p, Point v0, Point v1, Point v2){
	Point vec_a = v0 - p;
	Point vec_b = v1 - p;
	Point vec_c = v2 - p;

	double angle = solid_angle(vec_a, vec_b, vec_c);

	Point u = v1 - v0;
	Point v = v2 - v0;
	Point n = Point::cross(u, v);

	Point p_vec = p - v0;
	double dot = n * p_vec;

	int factor;
	if(dot > 0){
		factor = 1;
	} else {
		factor = -1;
	}

	return factor * angle;
}

double solid_angle(Point a, Point b, Point c){
	a.normalize();
	b.normalize();
	c.normalize();

	double numer = Point::cross(a, b) * c;
	double denom = 1 + (a * b) + (b * c) + (c * a);

	double angle = 2 * atan2(numer, denom);
	return fabs(angle);
}

}
}
