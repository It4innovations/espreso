
#include "meshpreprocessing.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"

#include <vector>

#include <fstream>
#include <sstream>
#include <iomanip>

#include <unistd.h>

namespace espreso {
namespace mesh {

bool edge_ray_intersect(Point p, Point v0, Point v1);

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

void computeVolumeIndices(ElementStore *elements, const NodeStore *nodes)
{
	profiler::syncstart("compute_volume_indices");

	// uniform grid
	esint grid_size = 21;
	std::vector<int> grid(grid_size * grid_size);
	double z = 0.0;
	/*Point grid_start = Point(-0.5, 0.5, z);
	Point grid_end = Point(0.5, -0.5, z);*/
//    Point grid_start = Point(-0.55, 0.55, z); // left + up
//    Point grid_end = Point(0.45, -0.45, z);
    Point grid_start = Point(-0.45, 0.45, z); // right + down
    Point grid_end = Point(0.55, -0.55, z);
	Point grid_offset = Point(((Point(grid_end.x, grid_start.y, z) - grid_start).length())/(grid_size - 1),
							((Point(grid_start.x, grid_end.y, z) - grid_start).length())/(grid_size - 1), z);
	printf("offset: %f %f %f\n", grid_offset.x, grid_offset.y, grid_offset.z);

	// elements cycle
	//usleep(20 * 1000000);
	int eindex = 0;
	for (auto e = elements->nodes->cbegin(); e != elements->nodes->cend(); ++e, ++eindex) {
		const Point &p_min_max = nodes->coordinates->datatarray()[*(e->begin())];
		Point el_min = Point(p_min_max.x, p_min_max.y, p_min_max.z);
		Point el_max = Point(p_min_max.x, p_min_max.y, p_min_max.z);

		for (auto n = e->begin(); n != e->end(); ++n) {
			Point &p = nodes->coordinates->datatarray()[*n];
			printf("%f %f %f\n", p.x, p.y, p.z);

			p.minmax(el_min, el_max); // min/max -> BB			
		}

		// grid points in BB
		int min_x_inx = (el_min.x - grid_start.x)/grid_offset.x + 1;
		int min_y_inx = (grid_start.y - el_max.y)/grid_offset.y + 1;
		int max_x_inx = (el_max.x - grid_start.x)/grid_offset.x;
		int max_y_inx = (grid_start.y - el_min.y)/grid_offset.y;

		for(int x = min_x_inx; x <= max_x_inx; x++){
			for(int y = min_y_inx; y <= max_y_inx; y++){
				// point in polygon
				Point p = Point(grid_start.x + x*grid_offset.x, grid_start.y - y*grid_offset.y, 0);
				int cn = 0;
				for (auto n = e->begin(); n != e->end() - 1; ++n) { // loop through edges
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
				grid[y*grid_size + x] = cn%2? eindex : 0;
				printf("cn %d\n", cn);
				printf("cell index %d\n", grid[y*grid_size + x]);

			}
		}

		// printf("min %f %f %f\n", el_min.x, el_min.y, el_min.z);
		// printf("max %f %f %f\n", el_max.x, el_max.y, el_max.z);
		printf("\n");
	}

	store(grid_size, grid);
	profiler::syncend("compute_volume_indices");
	eslog::checkpointln("MESH: VOLUME INDICES COMPUTED");
}

bool isPointInPolygon(ElementStore *elements, const NodeStore *nodes){

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
