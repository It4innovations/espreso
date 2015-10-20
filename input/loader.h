
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esmesh.h"
#include <string>

namespace esinput {

class ExternalLoader {

public:
	void load(mesh::Mesh &mesh)
	{
		points(mesh._coordinates);
		elements(mesh._elements);
		mesh._partPtrs.back() = mesh._elements.size();
		mesh.computeLocalIndices(0);
		boundaryConditions(mesh._coordinates);
		clusterBoundaries(mesh, mesh._clusterBoundaries);
	}

	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements) = 0;
	virtual void boundaryConditions(mesh::Coordinates &coordinates) = 0;
	virtual void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries) = 0;

	virtual ~ExternalLoader() {};
};

class InternalLoader {

public:
	void load(mesh::Mesh &mesh)
	{
		points(mesh._coordinates);

		elements(mesh._elements, mesh._partPtrs);
		for (size_t i = 0; i < mesh.parts(); i++) {
			mesh.computeLocalIndices(i);
		}

		fixPoints(mesh._fixPoints);
		size_t fSize = mesh._fixPoints.size() / mesh.parts();
		for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = 0; i < fSize; i++) {
				mesh._fixPoints[p * fSize + i] = mesh.coordinates().localIndex(mesh._fixPoints[p * fSize + i], p);
			}
		}

		mesh.computeBoundaries();
		corners(mesh._subdomainBoundaries);
		clusterBoundaries(mesh._clusterBoundaries);

		boundaryConditions(mesh._coordinates);
	}

	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void fixPoints(std::vector<eslocal> &fixPoints) = 0;
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
