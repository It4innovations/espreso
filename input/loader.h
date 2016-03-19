
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
		open();

		points(mesh._coordinates);
		elements(mesh._elements);
		faces(mesh._faces);
		boundaryConditions(mesh._coordinates);
		clusterBoundaries(mesh._clusterBoundaries, mesh._neighbours);

		close();

		partitiate(mesh._partPtrs);
		fixPoints(mesh._fixPoints);
		corners(mesh._subdomainBoundaries);
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
