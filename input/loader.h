
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esbasis.h"
#include "esmesh.h"
#include <string>

namespace espreso {

namespace input {

//class Loader

class ExternalLoader {

public:
	void load(Mesh &mesh)
	{
		open();
		points(mesh._coordinates);
		elements(mesh._elements);
		faces(mesh._faces);
		boundaryConditions(mesh._coordinates);
		clusterBoundaries(mesh, mesh._clusterBoundaries, mesh._neighbours);
		close();
		mesh.partitiate(config::mesh::subdomains);
		mesh.computeFixPoints(config::mesh::fixPoints);

		if (config::solver::FETI_METHOD == config::HYBRID_FETI) {
			mesh.computeCorners(
					config::mesh::corners,
					config::mesh::vertexCorners,
					config::mesh::edgeCorners,
					config::mesh::faceCorners,
					config::mesh::averageEdges,
					config::mesh::averageFaces);
		}
	}

protected:
	virtual void points(Coordinates &coordinates) = 0;
	virtual void elements(std::vector<Element*> &elements) = 0;
	virtual void faces(Faces &faces) = 0;
	virtual void boundaryConditions(Coordinates &coordinates) = 0;
	virtual void clusterBoundaries(Mesh &mesh, Boundaries &boundaries, std::vector<int> &neighbours) = 0;

	virtual void open() = 0;
	virtual void close() = 0;

	virtual ~ExternalLoader() {};
};

class APILoader {

public:
	void load(APIMesh &mesh)
	{
		points(mesh._coordinates);
		elements(mesh._elements);
		mesh.partitiate(config::mesh::subdomains);
		clusterBoundaries(mesh, mesh._clusterBoundaries, mesh._neighbours);

		mesh.computeCorners(
				config::mesh::corners,
				config::mesh::vertexCorners,
				config::mesh::edgeCorners,
				config::mesh::faceCorners,
				config::mesh::averageEdges,
				config::mesh::averageFaces);
	}

protected:
	virtual void points(Coordinates &coordinates) = 0;
	virtual void elements(std::vector<Element*> &elements) = 0;
	virtual void clusterBoundaries(Mesh &mesh, Boundaries &boundaries, std::vector<int> &neighbours) = 0;

	virtual ~APILoader() {};
};

class InternalLoader {

public:
	void load(Mesh &mesh)
	{
		points(mesh._coordinates);

		elements(mesh._elements, mesh._partPtrs);
		mesh.remapElementsToSubdomain();

		if (manualPartition()) {
			mesh.partitiate(mesh.parts());
			mesh.computeFixPoints(config::mesh::fixPoints);
		} else {
			fixPoints(mesh._fixPoints);
			for (size_t p = 0; p < mesh.parts(); p++) {
				for (size_t i = 0; i < mesh._fixPoints[p].size(); i++) {
					mesh._fixPoints[p][i] = mesh.coordinates().localIndex(mesh._fixPoints[p][i], p);
				}
				std::sort(mesh._fixPoints[p].begin(), mesh._fixPoints[p].end());
			}
			mesh.computeBoundaries();
		}
		boundaryConditions(mesh._coordinates);
		clusterBoundaries(mesh._clusterBoundaries, mesh._neighbours);

		if (manualPartition()) {
			mesh.computeCorners(
					config::mesh::corners,
					config::mesh::vertexCorners,
					config::mesh::edgeCorners,
					config::mesh::faceCorners,
					config::mesh::averageEdges,
					config::mesh::averageFaces);
		} else {
			corners(mesh._subdomainBoundaries);
			if (config::mesh::averageEdges || config::mesh::averageFaces) {
				mesh.computeCorners(0, true, false, false, config::mesh::averageEdges, config::mesh::averageFaces);
			}
		}
	}

protected:
	virtual bool manualPartition() = 0;

	virtual void points(Coordinates &coordinates) = 0;
	virtual void elements(std::vector<Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) = 0;
	virtual void boundaryConditions(Coordinates &coordinates) = 0;
	virtual void corners(Boundaries &boundaries) = 0;
	virtual void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours) = 0;

	virtual ~InternalLoader() {};
};

}
}



#endif /* INPUT_LOADER_H_ */
