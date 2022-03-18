
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

namespace espreso {
namespace mesh {

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

	esint voxels = 10;
	std::vector<int> grid(voxels * voxels);

	int eindex = 0;
	for (auto e = elements->nodes->cbegin(); e != elements->nodes->cend(); ++e, ++eindex) {
		for (auto n = e->begin(); n != e->end(); ++n) {
			const Point &p = nodes->coordinates->datatarray()[*n];
			printf("%f %f %f\n", p.x, p.y, p.z);
		}
		printf("\n");
	}

	store(voxels, grid);
	profiler::syncend("compute_volume_indices");
	eslog::checkpointln("MESH: VOLUME INDICES COMPUTED");
}

}
}
