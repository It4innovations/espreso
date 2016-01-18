
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esmesh.h"
#include <string>

namespace esinput {

class ExternalLoader {

public:
	void load(mesh::Mesh &mesh)
	{
		open();
		points(mesh._coordinates);
		elements(mesh._elements);
		mesh.partitiate(esconfig::mesh::subdomains);
		mesh.computeFixPoints(esconfig::mesh::fixPoints);
		mesh.computeBoundaries();
		boundaryConditions(mesh._coordinates);
		clusterBoundaries(mesh, mesh._clusterBoundaries);
		close();

		mesh.computeCorners(
				esconfig::mesh::corners,
				esconfig::mesh::vertexCorners,
				esconfig::mesh::edgeCorners,
				esconfig::mesh::faceCorners,
				esconfig::mesh::averaging);
	}

protected:
	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements) = 0;
	virtual void boundaryConditions(mesh::Coordinates &coordinates) = 0;
	virtual void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries) = 0;

	virtual void open() = 0;
	virtual void close() = 0;

	virtual ~ExternalLoader() {};
};

class InternalLoader {

public:
	void load(mesh::Mesh &mesh)
	{
		points(mesh._coordinates);

		elements(mesh._elements, mesh._partPtrs);
		mesh.remapElementsToSubdomain();

		if (manualPartition()) {
			mesh.partitiate(mesh.parts());
			mesh.computeFixPoints(esconfig::mesh::fixPoints);
		} else {
			fixPoints(mesh._fixPoints);
			for (size_t p = 0; p < mesh.parts(); p++) {
				for (size_t i = 0; i < mesh._fixPoints[p].size(); i++) {
					mesh._fixPoints[p][i] = mesh.coordinates().localIndex(mesh._fixPoints[p][i], p);
				}
				std::sort(mesh._fixPoints[p].begin(), mesh._fixPoints[p].end());
			}
		}
		mesh.computeBoundaries();
		boundaryConditions(mesh._coordinates);
		clusterBoundaries(mesh._clusterBoundaries);

		if (manualPartition()) {
			mesh.computeCorners(
					esconfig::mesh::corners,
					esconfig::mesh::vertexCorners,
					esconfig::mesh::edgeCorners,
					esconfig::mesh::faceCorners,
					esconfig::mesh::averaging);
		} else {
			corners(mesh._subdomainBoundaries);
		}
	}

protected:
	virtual bool manualPartition() = 0;

	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) = 0;
	virtual void boundaryConditions(mesh::Coordinates &coordinates) = 0;
	virtual void corners(mesh::Boundaries &boundaries) = 0;
	virtual void clusterBoundaries(mesh::Boundaries &boundaries) = 0;

	virtual ~InternalLoader() {};
};

template <class TLoader>
class Loader {

public:
	Loader(int argc, char** argv, int rank, int size): _loader(argc, argv, rank, size) { };

	void load(mesh::Mesh &mesh)
	{
		_loader.load(mesh);
	}

private:
	TLoader _loader;
};

}


#endif /* INPUT_LOADER_H_ */
