
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
		TimeEval measurement("Mesh loader"); measurement.totalTime.startWithBarrier();

		open();

		TimeEvent tPoints("coordinates"); tPoints.start();
		points(mesh._coordinates);
		tPoints.end(); measurement.addEvent(tPoints);

		TimeEvent tElements("elements"); tElements.start();
		elements(mesh._elements);
		tElements.end(); measurement.addEvent(tElements);

		TimeEvent tFaces("faces"); tFaces.start();
		faces(mesh._faces);
		tFaces.end(); measurement.addEvent(tFaces);

		TimeEvent tBoundaryConditions("boundary conditions"); tBoundaryConditions.start();
		boundaryConditions(mesh._coordinates);
		tBoundaryConditions.end(); measurement.addEvent(tBoundaryConditions);

		TimeEvent tClusterBoundaries("cluster boundaries"); tClusterBoundaries.start();
		clusterBoundaries(mesh._clusterBoundaries, mesh._neighbours);
		tClusterBoundaries.end(); measurement.addEvent(tClusterBoundaries);

		close();

		TimeEvent tPartition("partition"); tPartition.start();
		partitiate(mesh._partPtrs);
		tPartition.end(); measurement.addEvent(tPartition);

		TimeEvent tFixPoints("fix points"); tFixPoints.start();
		fixPoints(mesh._fixPoints);
		tFixPoints.end(); measurement.addEvent(tFixPoints);

		TimeEvent tCorners("corners"); tCorners.start();
		corners(mesh._subdomainBoundaries);
		tCorners.end(); measurement.addEvent(tCorners);

		measurement.totalTime.endWithBarrier(); measurement.printStatsMPI();
	}

protected:
	virtual void points(Coordinates &coordinates) = 0;
	virtual void elements(std::vector<Element*> &elements) = 0;
	virtual void faces(Faces &faces) { };
	virtual void boundaryConditions(Coordinates &coordinates) = 0;
	virtual void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours) = 0;

	virtual void open() {};
	virtual void close() {};

	virtual void partitiate(std::vector<eslocal> &parts)
	{
		mesh.partitiate(config::mesh::subdomains);
	}

	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints)
	{
		mesh.computeFixPoints(config::mesh::fixPoints);
	}

	virtual void corners(Boundaries &boundaries)
	{
		if (config::solver::FETI_METHOD == config::TOTAL_FETI) {
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
