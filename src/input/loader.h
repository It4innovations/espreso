
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esbasis.h"
#include "esmesh.h"
#include <string>

namespace espreso {
namespace input {

class Loader {

public:
	void fill()
	{
		auto sizeStats = [] (const std::vector<std::vector<eslocal> > &data) {
			std::vector<size_t> sizes(data.size());
			for (size_t p = 0; p < data.size(); p++) {
				sizes[p] = data[p].size();
			}
			return Info::averageValues(sizes);
		};
		auto intervalStats = [] (const std::vector<eslocal> &data) {
			std::vector<size_t> sizes(data.size() - 1);
			for (size_t p = 0; p < data.size() - 1; p++) {
				sizes[p] = data[p + 1] - data[p];
			}
			return Info::averageValues(sizes);
		};

		TimeEval measurement("Mesh loader"); measurement.totalTime.startWithBarrier();

		open();

		TimeEvent tPoints("coordinates"); tPoints.start();
		points(mesh._coordinates, mesh._DOFs);
		tPoints.end(); measurement.addEvent(tPoints);
		ESINFO(OVERVIEW) << "Coordinates loaded - total number of nodes: " << Info::sumValue(mesh.coordinates().clusterSize());

		TimeEvent tElements("elements"); tElements.start();
		elements(mesh._elements);
		materials(mesh._materials);
		tElements.end(); measurement.addEvent(tElements);
		ESINFO(OVERVIEW) << "Elements loaded - total number of elements: " << Info::sumValue(mesh.getElements().size());

		TimeEvent tFaces("faces"); tFaces.start();
		faces(mesh._faces);
		tFaces.end(); measurement.addEvent(tFaces);
		ESINFO(DETAILS) << "Faces loaded - total number of faces: " << Info::sumValue(mesh._faces.size());

		TimeEvent tSettings("settings"); tSettings.start();
		settings(mesh._evaluators, mesh._elements, mesh._coordinates);
		tSettings.end(); measurement.addEvent(tSettings);

		TimeEvent tClusterBoundaries("cluster boundaries"); tClusterBoundaries.start();
		clusterBoundaries(mesh._clusterBoundaries, mesh._neighbours);
		tClusterBoundaries.end(); measurement.addEvent(tClusterBoundaries);
		ESINFO(OVERVIEW) << "Neighbours loaded - number of neighbours for each cluster is " << Info::averageValue(mesh.neighbours().size());

		close();

		TimeEvent tPartition("partition"); tPartition.start();
		partitiate(mesh._partPtrs);
		tPartition.end(); measurement.addEvent(tPartition);
		ESINFO(OVERVIEW) << "Mesh partitioned into " << config::env::MPIsize << " * " << mesh.parts() << " = " << mesh.parts() * config::env::MPIsize
				<< " parts. There is " << intervalStats(mesh._partPtrs) << " elements in subdomain.";

		fixPoints(mesh._fixPoints);

		TimeEvent tCorners("corners"); tCorners.start();
		corners(mesh._subdomainBoundaries);
		tCorners.end(); measurement.addEvent(tCorners);

		measurement.totalTime.endWithBarrier(); measurement.printStatsMPI();
	}

protected:
	virtual void points(Coordinates &coordinates, size_t &DOFs) = 0;
	virtual void elements(std::vector<Element*> &elements) = 0;
	virtual void materials(std::vector<Material> &materials) = 0;
	virtual void faces(Faces &faces) { };
	virtual void settings(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, Coordinates &coordinates) {};
	virtual void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours) = 0;

	virtual void open() {};
	virtual void close() {};

	virtual void partitiate(std::vector<eslocal> &parts)
	{
		mesh.partitiate(config::mesh::subdomains);
	}

	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) {};

	virtual void corners(Boundaries &boundaries)
	{
		if (config::solver::FETI_METHOD == config::TOTAL_FETI || config::solver::B0_TYPE == config::B0Type::KERNELS) {
			return;
		}
		mesh.computeCorners(
				config::mesh::corners,
				config::mesh::vertexCorners,
				config::mesh::edgeCorners,
				config::mesh::faceCorners,
				config::mesh::averageEdges,
				config::mesh::averageFaces);
	}

	void remapElementsToSubdomains()
	{
		mesh.remapElementsToSubdomain();
	}

	void computeBoundaries()
	{
		mesh.computeBoundaries();
	}

	Loader(Mesh &mesh): mesh(mesh) {};
	virtual ~Loader() {};

protected:
	Mesh &mesh;
};

}
}



#endif /* INPUT_LOADER_H_ */
